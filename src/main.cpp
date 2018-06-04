#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define DBG_CSV 1
#define DBG_CONSOLE 0

#define NUM_NEXT_PTS_TO_SIM 50

// UNIT CONVERSIONS
#define KPH2MPS(x) (double)x / 3.6                             //
#define MPH2MPS(x) KPH2MPS((double)x*8./5.)                    //

#define TIME_STEP_sec (double)0.02                             // 20ms
#define MAX_TARGET_SPEED_MPS KPH2MPS((double)80.0 * (1-0.05))  // 80 kph * 3.6 mps/kph, 5% margin
#define LANE_WIDTH_M 4.0                                       // nominal lane width is 4.0m
#define GETLANE(d) ((int)d / LANE_WIDTH_M)                     // given d, return lane
#define MAX_ACCELERATION_MPSS (double)10.0 * (1.-0.1)          // 10 mpss given in requirements, 10% margin
#define MAX_JERK_MPSSS (double)50.0 * (1.-0.1)                 // 50 mpsss given in requirements, 10% margin

enum sensor_fusion_t { ID, PX, PY, VX, VY, S, D};
enum lane_t {left_lane, middle_lane, right_lane, num_lanes};

struct target_vehicle_struct {
  int ID;
  int lane;
  double deltaS = 9999; //initialize to arbitrary large number
  double longvel;
};

bool compare_target_vehicle_struct(const target_vehicle_struct &first, const target_vehicle_struct &second) {
  return (first.deltaS < second.deltaS);
}

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
#if(DBG_CSV)
  ofstream CSVFILE("output.csv");
  CSVFILE << "ego_px" << ", ";
  CSVFILE << "ego_py" << ", ";
  CSVFILE << "ego_s" << ", ";
  CSVFILE << "ego_d" << ", ";
  CSVFILE << "ego_heading_rad" << ", ";
  CSVFILE << "ego_longvel_mps" << ", ";
    for (int i = 0; i < NUM_NEXT_PTS_TO_SIM; i++) {
      CSVFILE << "unlim accel_demand" << i << ", ";
      CSVFILE << "future_ego_longvel_mps" << i << ", ";
      CSVFILE << "lim accel_demand" << i << ", ";
      CSVFILE << "future_x" << i << ", ";
//      CSVFILE << "future_y" << i << ", ";
    }
  CSVFILE << endl;
#endif

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &CSVFILE](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double ego_px = j[1]["x"];  // units: meters?
          double ego_py = j[1]["y"];  // units: meters?
          double ego_s = j[1]["s"];  // units: meters?
          double ego_d = j[1]["d"];  // units: meters?
          double ego_heading_rad = deg2rad(j[1]["yaw"]);   // units: radians?
          double ego_longvel_mps = MPH2MPS(j[1]["speed"]);  // units: m/s?
#if(DBG_CSV)
          CSVFILE << ego_px << ", ";
          CSVFILE << ego_py << ", ";
          CSVFILE << ego_s << ", ";
          CSVFILE << ego_d << ", ";
          CSVFILE << ego_heading_rad << ", ";
          CSVFILE << ego_longvel_mps << ", ";
#endif
          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];

          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          int ego_lane = GETLANE(ego_d);

          // GET NEAREST TARGET VEHICLE IN EACH LANE AND *AHEAD* OF EGO VEHICLE
          // TODO: CONSIDER FASTER TRAFFIC APPROACHING FROM BEHIND WHEN PLANNING LANE CHANGE
          vector < target_vehicle_struct> ClosestTrafficParticipant_in_lane(num_lanes);

          for (int i = 0; i < sensor_fusion.size(); i++) {
            target_vehicle_struct target_vehicle;

            target_vehicle.ID = sensor_fusion[i][ID];
            target_vehicle.lane = GETLANE(sensor_fusion[i][D]);
            target_vehicle.longvel = sqrt(pow((double)sensor_fusion[i][VX], 2) + pow((double)sensor_fusion[i][VY], 2));
            target_vehicle.deltaS = (double)sensor_fusion[i][S] - ego_s;

            if (target_vehicle.deltaS > 0) { // target in front of ego vehicle
              if (target_vehicle.deltaS < ClosestTrafficParticipant_in_lane[target_vehicle.lane].deltaS) {
                ClosestTrafficParticipant_in_lane[target_vehicle.lane] = target_vehicle;
              }
            }
          }
          //cout << "Lane " << ego_lane+1 << ":\t" << "ID: " << ClosestTrafficParticipant_in_lane[ego_lane].ID << "\tdeltaS: " << ClosestTrafficParticipant_in_lane[ego_lane].deltaS << endl;

          // Initialize goal state to stay in lane
          double goal_longvel_mps = MAX_TARGET_SPEED_MPS;
          int goal_lane = middle_lane;  

#if(0)
          // TODO: use state machine, build a path for each possible transition, calculate cost function, and then choose state with lowest cost
          if (left_lane == ego_lane) { //check straight and LCR
            //if ego slow
            if (   ego_longvel_mps < goal_longvel_mps \
                && ClosestTrafficParticipant_in_lane[middle_lane].longvel > ego_longvel_mps) { 
              goal_lane = middle_lane;
            }
          }
          else if (middle_lane == ego_lane) { // check straight, LCR, and LCL
            if (ego_longvel_mps < goal_longvel_mps \
              && ClosestTrafficParticipant_in_lane[right_lane].longvel > ego_longvel_mps) {
              goal_lane = right_lane;
            }
            if (ego_longvel_mps < goal_longvel_mps \
              && ClosestTrafficParticipant_in_lane[left_lane].longvel > ego_longvel_mps) {
              goal_lane = left_lane;
            }
          }
          else if (right_lane == ego_lane) { // check straight and LCL
            if (ego_longvel_mps < goal_longvel_mps \
              && ClosestTrafficParticipant_in_lane[middle_lane].longvel > ego_longvel_mps) {
              goal_lane = middle_lane;
            }
          }
          else { // this is an error state - just continue straight
          }
#endif
          //***********************************
          // Build the PLANNINGPATH 
          // Keep it simple with 4 points
          // Trust spline to do the smoothing
          //***********************************
          vector<double> planningpath_x_Fr0, planningpath_y_Fr0; // ego vehicle x/y-position in GLOBAL frame, Fr0
		      vector<double> planningpath_x_Fr1, planningpath_y_Fr1; // ego vehicle x/y-position in EGO_VEHICLE frame, Fr1

		      // POINT 1 - Include the vehicle's current position
          planningpath_x_Fr0.push_back(ego_px);
          planningpath_y_Fr0.push_back(ego_py);

          // POINT 2 - Project where ego vehicle will be in the next TIME_STEP_sec, if driving at 10. km/h
          planningpath_x_Fr0.push_back(ego_px + KPH2MPS(10.)*cos(ego_heading_rad)*TIME_STEP_sec);
          planningpath_y_Fr0.push_back(ego_py + KPH2MPS(10.)*sin(ego_heading_rad)*TIME_STEP_sec);

          // POINT 3 - Calculate the point 1 seconds ahead at max speed and in goal lane
          vector<double> goal_xy;
          goal_xy = getXY(ego_s + MAX_TARGET_SPEED_MPS*1.0, \
            (goal_lane * 4) + 2, \
            map_waypoints_s, \
            map_waypoints_x, \
            map_waypoints_y);
          planningpath_x_Fr0.push_back(goal_xy[0]);
          planningpath_y_Fr0.push_back(goal_xy[1]);

          // POINT 4 - Calculate the point 2.5 seconds ahead at max speed and in goal lane
          //vector<double> goal_xy;
          goal_xy = getXY(ego_s + MAX_TARGET_SPEED_MPS*2.5, \
                     (goal_lane*4) + 2, \
                     map_waypoints_s, \
                     map_waypoints_x, \
                     map_waypoints_y);
          planningpath_x_Fr0.push_back(goal_xy[0]);
          planningpath_y_Fr0.push_back(goal_xy[1]);

          // transform from global reference frame to ego vehicle reference frame
		      // https://en.wikipedia.org/wiki/Rotation_of_axes
          for (int i = 0; i < planningpath_x_Fr0.size(); i++) {
            // shift
			      double shift_x = planningpath_x_Fr0[i] - ego_px;
			      double shift_y = planningpath_y_Fr0[i] - ego_py;
            // rotate
            planningpath_x_Fr1.push_back( shift_x * cos(ego_heading_rad) + shift_y * sin(ego_heading_rad));
            planningpath_y_Fr1.push_back(-shift_x * sin(ego_heading_rad) + shift_y * cos(ego_heading_rad));
          }
#if(0)
          cout << planningpath_x_Fr1[0] << "\t" << planningpath_y_Fr1[0] << "\t" \
               << planningpath_x_Fr1[1] << "\t" << planningpath_y_Fr1[1] << "\t" \
               << planningpath_x_Fr1[2] << "\t" << planningpath_y_Fr1[2] << "\t" \
               << planningpath_x_Fr1[3] << "\t" << planningpath_y_Fr1[3] << "\t" \
               << endl;
#endif
          //*********************************
          // Build the planning_path_spline
          //*********************************
          tk::spline planning_path_spline;
          planning_path_spline.set_points(planningpath_x_Fr1, planningpath_y_Fr1);

          // ********************
          // Build path request
          // ********************
          vector<double> next_x_vals, next_y_vals; // ego vehicle desired x/y-position in GLOBAL frame, Fr0, to send to motion control / simulator

#if(1)
		      // Load up the planning path with the unused points from the previous path
		      for (int i = 0; i < previous_path_x.size(); i++) {
			      planningpath_x_Fr0.push_back(previous_path_x[i]);
			      planningpath_y_Fr0.push_back(previous_path_y[i]);
		      }
#endif

          double future_ego_longvel_mps = ego_longvel_mps;
          double accel_demand=0., accel_demand_z1=0.;
          double future_x_z1 = 0.;
#if(1)
          cout << previous_path_x.size() << endl;
          for (int i = 0; i < (NUM_NEXT_PTS_TO_SIM - previous_path_x.size()); i++) {
#else
          for (int i = 0; i < NUM_NEXT_PTS_TO_SIM; i++) {
#endif
            if (future_ego_longvel_mps < goal_longvel_mps) { // need to accelerate
			        accel_demand = MAX_ACCELERATION_MPSS;
            } else if (future_ego_longvel_mps > goal_longvel_mps) { // need to decelerate
			        accel_demand = -MAX_ACCELERATION_MPSS;
            } else { 
              // maintain speed
              accel_demand = 0.;
            }
#if(DBG_CONSOLE)
			cout << "accel_demand: " << accel_demand << "\t";
#endif
#if(DBG_CSV)
      CSVFILE << accel_demand << ", ";
#endif
            //limit jerk            
            accel_demand = TIME_STEP_sec * max(-MAX_JERK_MPSSS, min(MAX_JERK_MPSSS, (accel_demand - accel_demand_z1) / TIME_STEP_sec));
            future_ego_longvel_mps += accel_demand * TIME_STEP_sec;
            accel_demand_z1 = accel_demand;
#if(DBG_CONSOLE)
			cout << "lim. accel. demand: " << accel_demand << "\t";
#endif
#if(DBG_CSV)
      CSVFILE << future_ego_longvel_mps << ", ";
      CSVFILE << accel_demand << ", ";
#endif
            double future_x, future_y; // future ego xy-position in EGO_VEHICLE frame, Fr1
            future_x = future_x_z1 + future_ego_longvel_mps * TIME_STEP_sec;
            future_x_z1 = future_x;
            future_y = planning_path_spline(future_x);
#if(DBG_CONSOLE)
			cout << "future_x: " << future_x << "\t";
			cout << "future_y: " << future_y << "\t";
#endif
#if(DBG_CSV)
      CSVFILE << future_x << ", ";
//      CSVFILE << future_y << ", ";
#endif
            // swing back to GLOBAL reference frame, Fr0
            // rotate
            future_x = future_x * cos(ego_heading_rad) - future_y * sin(ego_heading_rad);
            future_y = future_x * sin(ego_heading_rad) + future_y * cos(ego_heading_rad);

            // shift
            future_x += ego_px;
            future_y += ego_py;

            next_x_vals.push_back(future_x);
            next_y_vals.push_back(future_y);
          }
#if(DBG_CONSOLE)
		  cout << endl << endl;
#endif
#if (DBG_CSV)
      CSVFILE << endl;
#endif
          //END
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
