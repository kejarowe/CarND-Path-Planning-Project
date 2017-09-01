#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <map>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define LANE_WIDTH 4.0

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

struct map_data {
  std::vector<double> * x;
  std::vector<double> * y;
  std::vector<double> * s;
};

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
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
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

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
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

inline void next_xy(double target_speed, double target_d, double &pos_x, double &pos_y, double &pos_theta, map_data &md, double dt = 0.02)
{
  static double current_speed = 0;
  static double current_accel = 0;
  static double current_theta = 0;
  static double max_d_theta = pi(); 
  static double allowed_lon_accel = 5;
  static double allowed_jerk = 5;
  const  double allowed_lat_accel = 5;
  const  double projection_multiplier = 5; //10 //5
  const  double projection_offset = 5; //15

  //limit centripital acceleration
  if (current_speed > 0)
    max_d_theta = allowed_lat_accel * dt / current_speed;

  double speed_delta = target_speed - current_speed;
  allowed_lon_accel = std::copysign(allowed_lon_accel,speed_delta);
  allowed_jerk = std::copysign(allowed_jerk,speed_delta);

  double current_dist = current_speed * dt;
  double look_ahead_distance = current_dist * projection_multiplier + projection_offset; 
  
  //set jerk
  double desired_accel = 0; 
  if (std::abs(1.0 / dt * speed_delta) <= allowed_lon_accel ){
    //we can reach our desired velocity in the next time step
    desired_accel = 1.0 / dt * speed_delta;
  }
  else if (std::abs(speed_delta) <= 0.5 * std::pow(current_accel,2) / allowed_jerk) {
    //max jerk towards 0 so that we reach our speed our accel will be 0
    desired_accel = 0;
  }
  else {
    //max acceleration in the direction of desired speed
    desired_accel = std::copysign(allowed_lon_accel,speed_delta);
  }

  double delta_accel = desired_accel - current_accel;
  double jerk_command = std::copysign(std::min(allowed_jerk,std::abs(delta_accel)),delta_accel);


  current_speed += current_accel * dt;
  current_accel += jerk_command * dt;
  
  vector<double> sd_vec = getFrenet(pos_x, pos_y, pos_theta, *md.x, *md.y);
  vector<double> next_xy = getXY(sd_vec[0] + look_ahead_distance, target_d, *md.s, *md.x, *md.y);

  current_theta = std::atan2(next_xy[1] - pos_y, next_xy[0] - pos_x);
  double d_theta = current_theta - pos_theta;
  //fix angle wrap around
  while (d_theta > M_PI)
    d_theta -= 2*M_PI;
  while(d_theta < -M_PI)
    d_theta += 2*M_PI;
  pos_theta += std::abs(d_theta) > max_d_theta ? std::copysign(max_d_theta,d_theta) : d_theta;
  pos_x += current_dist * cos(pos_theta);
  pos_y += current_dist * sin(pos_theta);
}

int lane_from_d(double d)
{
  return std::round( d/LANE_WIDTH - 0.5 );
}

double magnitude(double a, double b)
{
  return std::sqrt(std::pow(a,2)+std::pow(b,2));
}

bool cipv_and_lane_speeds(double ego_s, double ego_d, std::vector<std::vector<double>> &car_list, double &cipv_ds, double &cipv_speed, double &left_lane_speed, double &right_lane_speed)
{
  const double safe_traffic_distance = 30;

  int ego_lane = lane_from_d(ego_d);
  
  left_lane_speed = ego_lane == 0? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
  right_lane_speed = ego_lane == 2? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
  
  
  cipv_ds = std::numeric_limits<double>::infinity();
  bool result = false;

  int num_right_lane_cars = 0;
  int num_left_lane_cars = 0;
  
  for (auto &car : car_list) {

    int car_lane = lane_from_d(car[6]);
    double ds = car[5] - ego_s;

    if (car_lane == ego_lane) {
      //car is in our lane, find cipv distance and speed
      if (ds > 0) {
        //in front of ego car
        if (ds < cipv_ds) {
          //this car is closer than any car we have seen yet
          result = true;
          cipv_ds = ds;
          cipv_speed = magnitude(car[3],car[4]);
        }
      }
    }
    else if (std::abs(car_lane - ego_lane) == 1){
      /*car is in adjacent lane
        -if its too close, set lane speed to -infinity
        -if its far enough ahead, set lane speed to min(car speed, lane_speed)
      */
      double &ref_lane_speed = car_lane > ego_lane ? right_lane_speed : left_lane_speed;
      car_lane > ego_lane ? num_right_lane_cars++ : num_left_lane_cars++;

      if (std::abs(ds) < safe_traffic_distance) {
        ref_lane_speed = -std::numeric_limits<double>::infinity();
      }
      else if (ds > 0) {
        ref_lane_speed = std::min(ref_lane_speed, magnitude(car[3],car[4]));
      }
    }
  }
  //printf("left speed: %lf, right_speed: %lf\n",left_lane_speed,right_lane_speed);
  //printf("%d left cars, %d right cars\n",num_left_lane_cars,num_right_lane_cars); 
  return result;
}

void traj_goals(double ego_d, double target_vehicle_speed, double target_vehicle_ds, double left_lane_speed, double right_lane_speed, double &desired_ego_speed, double &desired_ego_d)
{
  const double psuedo_target_headway_threshold = 2.5; //seconds, 5 was too long
  const double set_speed = 20.1; //45mph in mps
  const double lane_change_threshold_factor = 0.9;
  static int desired_lane = 1; //leftmost lane: 0, middle_lane: 1, right_lane: 2

  int ego_lane = lane_from_d(ego_d);
  bool changing_lanes = ego_lane != desired_lane;
  double psuedo_target_headway = target_vehicle_ds / set_speed;

  if (psuedo_target_headway < psuedo_target_headway_threshold) {
    //match target speed since its close to us
    desired_ego_speed = std::min(target_vehicle_speed,set_speed);
  }
  else {
    //no nearby targets in our lane
    desired_ego_speed = set_speed;
  }

  /*
    - if desired ego speed < 0.9 set_speed try to changes lanes, unless not in desired lane
  */

  if (!changing_lanes && desired_ego_speed < lane_change_threshold_factor * set_speed) {
    //our speed is too slow and we are not changing lanes, look for available lane

    if (left_lane_speed > right_lane_speed && left_lane_speed > desired_ego_speed) {
      //left lane is faster than our lane and right lane
      desired_lane = ego_lane - 1;
      printf("lane speed is %lf, left lane speed is: %lf, changing\n",desired_ego_speed,left_lane_speed);
    }
    else if (right_lane_speed >= left_lane_speed && right_lane_speed > desired_ego_speed){
      //right lane is faster than our lane and left lane
      desired_lane = ego_lane + 1;
      printf("lane speed is %lf, right lane speed is: %lf, changing\n",desired_ego_speed,right_lane_speed);
    }
  }

  //convert lane to d value
  desired_ego_d = (desired_lane + 0.5) * LANE_WIDTH;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  map_data md;
  md.x = &map_waypoints_x;
  md.y = &map_waypoints_y;
  md.s = &map_waypoints_s;

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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&md](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            static double last_speed = car_speed;
            std::chrono::time_point<std::chrono::system_clock> current_time = std::chrono::system_clock::now();
            static std::chrono::time_point<std::chrono::system_clock> last_time = current_time;
            double dt = (current_time - last_time).count();
            double car_accel = dt > 0 ? (car_speed - last_speed) / dt : 0;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            std::vector<std::vector<double>> car_list = sensor_fusion;

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            /*
              Trajectory Planning Strategy
              -find cipv by parsing fused object list
              -send cipv to decision maker which will determine desired speed and d
              -send desired speed and d to next_xy
              -next_xy calculates point to send to dynamic controller
            */

            //trajectory generation vars
            double pos_x, pos_y, pos_theta, pos_s, pos_d,
              pos_speed, pos_accel, desired_d, desired_speed;
            std::vector<double> spline_input_x;
            std::vector<double> spline_input_y;
            
            //copy old trajectory
            int path_size = previous_path_x.size();
            for(int i = 0; i < path_size; i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }
            
            //get last point on old trajectory
            if(path_size == 0)
            {
              pos_x = car_x;
              pos_y = car_y;
              pos_s = car_s;
              pos_d = car_d;
              pos_theta = 0;
              pos_speed = car_speed;
              pos_accel = car_accel;
            }
            else
            {
              pos_x = next_x_vals[path_size-1];
              pos_y = next_y_vals[path_size-1];
              if (path_size > 1)
                pos_theta = std::atan2(pos_y - next_y_vals[path_size-2], pos_x - next_x_vals[path_size-2]);
              else
                pos_theta = 0;
              vector<double> sd_vec = getFrenet(pos_x, pos_y, pos_theta, map_waypoints_x, map_waypoints_y);
              pos_s = sd_vec[0];
              pos_d = sd_vec[1];
              pos_speed = car_speed;
              pos_accel = car_accel;
            }

            //generate new trajectory points

            //find cipv and lane speeds
            double target_vehicle_ds, target_vehicle_speed, left_lane_speed, right_lane_speed;
            bool success = cipv_and_lane_speeds(pos_s, pos_d, car_list, target_vehicle_ds, target_vehicle_speed, left_lane_speed, right_lane_speed); 
        
            //run fsm
            traj_goals(pos_d, target_vehicle_speed, target_vehicle_ds, left_lane_speed, right_lane_speed, desired_speed, desired_d);

            double next_x, next_y;
            next_x = pos_x;
            next_y = pos_y;

            //generate next points
            for(int i = 0; i < 50-path_size; i++)
            { 
              next_xy(desired_speed, desired_d, pos_x, pos_y, pos_theta, md);
              //convert pos
              next_x_vals.push_back(pos_x);
              next_y_vals.push_back(pos_y);
            }

            last_speed = car_speed;
            last_time = current_time;
            //end path generation

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
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
















































































