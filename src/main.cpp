#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/QR>
#include "json.hpp"
#include "jmt.h"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

//Convert MPH to KPH
double mph2kph(double mph)
{
  return mph*1.60934;
}
double kph2mph(double kph)
{
  return kph/1.60934;
}

double kph2mps(double kph)
{
  return kph*0.27777778;
}

int num_samples=100;
int lane=1;
double ref_vel = kph2mps(mph2kph(40));


std::map<int, double[4]> cars;

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

void Planner(std::map<int, double[4]> cars, double car_s, int &lane, double &speed)
{
  double lane_speed[3]={-1,-1,-1};
  double lane_dist[3]={999,999,999};
  double lane_cost[3]={0,0,0};

  for (auto it=cars.begin();it!=cars.end();it++)
  {
    double s = it->second[0] - car_s;
    int l = it->second[3];
    double spd = it->second[2];

    printf("l %d  s %f  spd %f\n", l,s,spd);
    if(fabs(s)>30)
    {
      lane_speed[l]=max(lane_speed[l],spd);
    }

    lane_dist[l]=min(lane_dist[l],fabs(s));

    //slow down
    if(l==lane && s>0 && s<50)
    {
      speed = spd-kph2mps(5);
      // return;
    }
  }
  printf("Lane speed: %f %f %f\t", lane_speed[0],lane_speed[1],lane_speed[2]);
  printf("dist: %f %f %f\t", lane_dist[0],lane_dist[1],lane_dist[2]);


  int max_spd_lane=0;
  int min_spd_lane=0;
  int max_dist_lane=0;
  int min_dist_lane=0;
  for(int i=0;i<3;++i)
  {
    //Lane change cost
    if (i!=lane)
    {
      lane_cost[i]+=150;

      if(abs(lane-i)>1)
        lane_cost[i]+=1000;

    }

    if(lane_dist[i]<10)
    {
      lane_cost[i]+=1000;
    }
    if(lane_dist[i]<30)
    {
      lane_cost[i]+=500;
    }
    if(lane_dist[i]<50)
    {
      lane_cost[i]+=100;
    }

    lane_cost[i]+=50-lane_speed[i];


    if(lane_speed[i]>lane_speed[max_spd_lane])
    {
      max_spd_lane=i;
    }
    if(lane_speed[i]<lane_speed[min_spd_lane])
    {
      min_spd_lane=i;
    }

    if(lane_dist[i]>lane_dist[max_dist_lane])
    {
      max_dist_lane=i;
    }
    if(lane_dist[i]<lane_dist[min_dist_lane])
    {
      min_dist_lane=i;
    }
  }
  lane_cost[min_spd_lane]+=100;
  lane_cost[min_dist_lane]+=100;

  printf("cost: %f %f %f\n", lane_cost[0],lane_cost[1],lane_cost[2]);

  int lowest_cost_lane=lane;

  for(int i=0;i<3;++i)
  {
    if (lane_cost[i]<lane_cost[lowest_cost_lane])
      lowest_cost_lane = i;
  }
  lane = lowest_cost_lane;
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

  double car_yaw = abs(theta-heading);

  if(car_yaw > pi()/4)
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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
            double car_yaw = deg2rad(j[1]["yaw"]);
            double car_speed = j[1]["speed"];

            // Previous path data given to the Planner
            auto previous_path_x = j[1]["previous_path_x"];
            auto previous_path_y = j[1]["previous_path_y"];
            // Previous path's end s and d values
            double end_path_s = j[1]["end_path_s"];
            double end_path_d = j[1]["end_path_d"];

            // Sensor Fusion Data, a list of all other cars on the same side of the road.
            auto sensor_fusion = j[1]["sensor_fusion"];

            json msgJson;



            for (uint i=0;i<sensor_fusion.size();++i)
            {
              // std::cout << sensor_fusion[i] <<'\n';
              int id = sensor_fusion[i][0];

              cars[id][0] = sensor_fusion[i][5]; //s
              cars[id][1] = sensor_fusion[i][6]; //d
              cars[id][2] = sqrt(double(sensor_fusion[i][3])*double(sensor_fusion[i][3]) + double(sensor_fusion[i][4])*double(sensor_fusion[i][4])); //v
              cars[id][3] = int((double(sensor_fusion[i][6])-2.0)/4.0+0.5);

              // printf("Lane: %f Dist: %f speed: %f\n", cars[id][3],cars[id][0]-car_s,cars[id][2]);
            }


            vector<double> next_x_vals;
            vector<double> next_y_vals;

            vector<double> ptsx;
            vector<double> ptsy;

            double curr_vel = kph2mps(mph2kph(car_speed));

            double last_end_x;
            double last_end_y;
            int path_size = previous_path_x.size();

            //Add all previous points
            for(int i = 0; i < path_size; i++)
            {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }

            Planner(cars,car_s,lane,ref_vel);

            if(path_size < 2)
            {
                last_end_x = car_x;
                last_end_y = car_y;
                end_path_s = car_s;

                double last_end_x2 = car_x - cos(car_yaw);
                double last_end_y2 = car_y - sin(car_yaw);

                ptsx.push_back(last_end_x2);
                ptsx.push_back(last_end_x);

                ptsy.push_back(last_end_y2);
                ptsy.push_back(last_end_y);
            }
            else
            {
                last_end_x = previous_path_x[path_size-1];
                last_end_y = previous_path_y[path_size-1];

                double last_end_x2 = previous_path_x[path_size-2];
                double last_end_y2 = previous_path_y[path_size-2];
                car_yaw = atan2(last_end_y-last_end_y2,last_end_x-last_end_x2);

                ptsx.push_back(last_end_x2);
                ptsx.push_back(last_end_x);

                ptsy.push_back(last_end_y2);
                ptsy.push_back(last_end_y);
            }
            double dist_inc = 0.3;

            vector< double> start ={end_path_s,curr_vel,0};
            double T=(1+num_samples-path_size)*0.02;
            vector <double> end={end_path_s+T*ref_vel,ref_vel,0};
            std::vector<double> traj = JMT(start, end, T);

            // printf("%f,%f,%f => %f,%f,%f\t", start[0], start[1], start[2],end[0],end[1],end[2]);
            // printf("curr_vel: %f  ", curr_vel);
            // printf("Traj: ");
            // for(int i=0;i<traj.size();++i)
            // {
            //   printf("%f ", i,traj[i]);
            // }
            // printf("\n");




            // printf("Car %f  End %f\n", car_s, end_path_s);

            std::vector<double> next_wp0 = getXY(end_path_s+30,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            std::vector<double> next_wp1 = getXY(end_path_s+60,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            std::vector<double> next_wp2 = getXY(end_path_s+90,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            std::vector<double> next_wp3 = getXY(end_path_s+120,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);
            ptsx.push_back(next_wp3[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);
            ptsy.push_back(next_wp3[1]);

            // for(int i=0;i<ptsx.size();++i)
            // {
            //     printf("pts: %f, %f\n", ptsx[i], ptsy[i]);
            // }

            printf("-------------------\n");
            printf("%f %f %f\n", last_end_x,last_end_y,car_yaw);
            printf("-------------------\n");
            for(int i=0;i<ptsx.size();++i)
            {
                // printf("%f\t%f\n", ptsx[i], ptsy[i]);
                double shift_x = ptsx[i]-last_end_x;
                double shift_y = ptsy[i]-last_end_y;

                ptsx[i] = (shift_x * cos(0) - shift_y*sin(0));
                ptsy[i] = (shift_y * sin(0) + shift_y*cos(0));

                // ptsx[i] = (shift_x * cos(0-car_yaw) - shift_y*sin(0-car_yaw));
                // ptsy[i] = (shift_y * sin(0-car_yaw) + shift_y*cos(0-car_yaw));
                // printf("\t\t%f\t%f\n", ptsx[i], ptsy[i]);

            }

            // printf("------------\n");

            // Create spline
            tk::spline s;

            // for(int i=0;i<ptsx.size();++i)
            // {
            //     printf("pts: %f, %f\n", ptsx[i], ptsy[i]);
            // }
            //
            // printf("===========================================\n");


            s.set_points(ptsx,ptsy);
            //
            double tmp_last_s=0;
            //
            double target_x = 10;
            double target_y = s(target_x);
            double target_dist = sqrt(target_x*target_x + target_y*target_y);
            double x_add_on = 0;

            // for(int i = 0; i < num_samples-path_size; i++)
            // {
            //     // ========================== Good working JMT ==============================
            //     double t=(i)*0.02;
            //     double t2=t*t;
            //     double t3=t2*t;
            //     double t4=t3*t;
            //     double t5=t4*t;
            //     double next_s = traj[0] + traj[1]*t + traj[2]*t2 + traj[3]*t3 + traj[4]*t4 + traj[5]*t5;
            //     // printf("next_s %f  last + %f\n", next_s, next_s-tmp_last_s);
            //     tmp_last_s = next_s;
            //     // double next_s = end_path_s + (i+1)*dist_inc;
            //     double next_d = 2+4*lane;
            //     // printf("map_waypoints size: %d %d %d\n", map_waypoints_s.size(),map_waypoints_x.size(),map_waypoints_y.size());
            //     vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            //
            //     // printf("\t\t%f\t%f\n", xy[0],xy[1]);
            //
            //     next_y_vals.push_back(xy[1]);
            //     next_x_vals.push_back(xy[0]);
            //   }
              for(int i = 0; i < num_samples-path_size; i++)
              {
                // ======================================= Spline =======================================
                double N = target_dist/(.02*ref_vel);
                double x_point = x_add_on + target_x/N;
                double y_point = s(x_point);

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                x_point = x_ref*cos(0) - y_ref*sin(0);
                y_point = x_ref*sin(0) + y_ref*cos(0);

                x_point += last_end_x;
                y_point += last_end_y;

                // printf("spline %f, %f\n", x_point, y_point);
                // printf("\t\t\t\t%f\t%f\n", x_point,y_point);

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }


            msgJson["next_x"] = next_x_vals;
            msgJson["next_y"] = next_y_vals;

            auto msg = "42[\"control\","+ msgJson.dump()+"]";

            this_thread::sleep_for(chrono::milliseconds(100));
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
