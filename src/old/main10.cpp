#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

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

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Evaluate a polynomial derivative
double polyderv(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += coeffs[i]*i*CppAD::pow(x, i-1);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
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

// Transform from Frenet s,d coordinates to Cartesian x,y using spline
vector<double> getXY(double s, double d, tk::spline wpx, tk::spline wpy, tk::spline wpdx, tk::spline wpdy) {
  double x = wpx(s) + d*wpdx(s);
  double y = wpy(s) + d*wpdy(s);
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

  // for our spline fit
  vector<double> vx;
  vector<double> vy;
  vector<double> vd;
  vector<double> xyd;
  tk::spline wpx;
  tk::spline wpy;
  tk::spline wpdx;
  tk::spline wpdy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  double nextd = 6.;
  double dist_inc = 0.442;
  int timestep = 0;
  bool lanechange = false;

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

  // add extra waypoint to loop back to the beginning with overlap...
  map_waypoints_x.push_back(map_waypoints_x[0]);
  map_waypoints_y.push_back(map_waypoints_y[0]);
  map_waypoints_s.push_back(max_s);
  map_waypoints_dx.push_back(map_waypoints_dx[0]);
  map_waypoints_dy.push_back(map_waypoints_dy[0]);
  map_waypoints_x.push_back(map_waypoints_x[1]);
  map_waypoints_y.push_back(map_waypoints_y[1]);
  map_waypoints_s.push_back(max_s+map_waypoints_s[1]);
  map_waypoints_dx.push_back(map_waypoints_dx[1]);
  map_waypoints_dy.push_back(map_waypoints_dy[1]);

  // set up the spline
  wpx.set_points(map_waypoints_s,map_waypoints_x);
  wpy.set_points(map_waypoints_s,map_waypoints_y);
  wpdx.set_points(map_waypoints_s,map_waypoints_dx);
  wpdy.set_points(map_waypoints_s,map_waypoints_dy);

  // set up logging
  string log_file = "../data/logger.csv";
  ofstream out_log(log_file.c_str(), ofstream::out);
  out_log << "s,d,vd,x,y,xyd" << endl;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&wpx,&wpy,&wpdx,&wpdy,&vx,&vy,&vd,&xyd,&nextd,&dist_inc,&out_log,&timestep,&lanechange](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int path_size = previous_path_x.size();
          int num_points = 100;
          double angle = deg2rad(car_yaw);
          double pos_x = car_x;
          double pos_y = car_y;
          double max_s = 6945.554;
          vector<double> t;
          vector<double> t2;
          vector<double> inc;
          vector<double> nd;
          vector<double> frenet;
          tk::spline smooth_speed;
          tk::spline smooth_local;
          tk::spline smooth_lanechange;

          // at beginning - no paths
          if (path_size == 0)
          {
            t.push_back(0.);
            t.push_back(double(num_points*5));
            inc.push_back(0.025);
            inc.push_back(dist_inc);
            smooth_speed.set_points(t,inc);

            nd.push_back(nextd);
            nd.push_back(nextd);
            smooth_lanechange.set_points(t,nd);

            frenet = getFrenet(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
            double nexts = frenet[0];
            for (int i = 0; i<num_points; i++)
            {
              nexts += smooth_speed(double(i));
              vx.push_back(nexts);
              // vy.push_back(frenet[1]);
              vy.push_back(nextd);
              if (i > 0)
                vd.push_back(distance(vx[i-1], vy[i-1], vx[i], vy[i]));
              else
                vd.push_back(0.);
              vector<double> xy2 = getXY(nexts, nextd, wpx, wpy, wpdx, wpdy);
              next_x_vals.push_back(xy2[0]);
              next_y_vals.push_back(xy2[1]);
            }
          }

          // we are already moving...
          else
          {
            pos_x = previous_path_x[path_size-1];
            pos_y = previous_path_y[path_size-1];

            double pos_x2 = previous_path_x[path_size-2];
            double pos_y2 = previous_path_y[path_size-2];
            angle = atan2(pos_y-pos_y2,pos_x-pos_x2);

            // std::cout << "pathsize: " << path_size << " , " << vx.size() << std::endl;
            for (int i = 0; i < path_size; i++)
            {
              vector<double> xy = getXY(vx[i], vy[i], wpx, wpy, wpdx, wpdy);
              // cout << i << " s,d: " <<  vx[i] << "," << vy[i] << " vd: " << vd[i] << " x,y: " << xy[0] << "," << xy[1] << " xyd: " << xyd[i] << std::endl;
            }

            for (int i = 0; i < (num_points-path_size); i++)
            {
              vector<double> xy = getXY(vx[i], vy[i], wpx, wpy, wpdx, wpdy);
              out_log << vx[i] << "," << vy[i] << "," << vd[i] << "," << xy[0] << "," << xy[1] << "," << xyd[i] << std::endl;
            }
            vx.erase(vx.begin(),vx.begin()+(num_points-path_size));
            vy.erase(vy.begin(),vy.begin()+(num_points-path_size));
            vd.erase(vd.begin(),vd.begin()+(num_points-path_size));
            xyd.erase(xyd.begin(),xyd.begin()+(num_points-path_size));

            frenet.push_back(vx[vx.size()-1]);
            frenet.push_back(vy[vy.size()-1]);
            if (lanechange && abs(vy[0]-nextd) < 0.1)
              lanechange = false;

            for(int i = 0; i < path_size; i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            t.push_back(0.);
            t.push_back(double(num_points*5));
            if (vd[0] <= dist_inc)
              inc.push_back(vd[0]);
            else
              inc.push_back(dist_inc);
            inc.push_back(dist_inc);
            smooth_speed.set_points(t,inc);

            t2.push_back(double(0));
            t2.push_back(double(path_size-1));
            t2.push_back(double(num_points*1.5));
            nd.push_back(vy[0]);
            nd.push_back(vy[path_size-1]);
            nd.push_back(nextd);
            smooth_lanechange.set_points(t2,nd);

            // std::cout << "x,y: " << pos_x << "," << pos_y << " frenet: " << frenet[0] << "," << frenet[1] << std::endl;

            double nexts = frenet[0];
            for (int i = path_size; i<num_points; i++)
            {
              nexts += smooth_speed(double(i));
              if (nexts > max_s)
                nexts -= max_s;
              vx.push_back(nexts);
              // vy.push_back(frenet[1]);
              vy.push_back(smooth_lanechange(double(i)));
              if (i > 0)
                vd.push_back(distance(vx[i-1], vy[i-1], vx[i], vy[i]));
              else
                vd.push_back(0.);
              vector<double> xy2 = getXY(vx[i], vy[i], wpx, wpy, wpdx, wpdy);
              next_x_vals.push_back(xy2[0]);
              next_y_vals.push_back(xy2[1]);
            }
          }

          if (timestep > 50 && not lanechange)
          {
            frenet = getFrenet(car_x, car_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);
            cout << "frnet: " << frenet[1] << " current position: " << frenet[0] << "," << frenet[1] << endl;
            cout << "nextd: " << nextd << " current position: " << vx[0] << "," << vy[0] << endl;
            cout << "currd: " << vy[0] << " " << path_size << " position: " << vx[path_size-1] << "," << vy[path_size-10] << endl;
            for (int k = 0; k<sensor_fusion.size(); k++)
            {
              auto vid = sensor_fusion[k];
              double vidvx = vid[3];
              double vidvy = vid[4];
              double vids = vid[5];
              double vidd = vid[6];
              cout << "vid: " << vid[0] << " current position: " << vids << "," << vidd << " spotted at: ";
              for (int i = 0; i<vx.size(); i++)
              {
                if (abs(vx[i]-vids) < 5.0 && abs(nextd-vidd) < 1.0)
                {
                  cout << i << " ";
                  if ( abs(nextd-vy[i]) < 1.0 and not lanechange )
                  {
                    if ( nextd < 7. && i < 80)
                      nextd += 4.;
                    else if (nextd > 5.)
                      nextd -= 4.;
                    lanechange = true;
                  }
                  /*
                  if ( abs(nextd-vy[i]) < 1.0)
                  {
                    dist_inc = sqrt(vidvx*vidvx+vidvy*vidvy)/50.;
                  }
                  */
                }
              }
              cout << endl;
            }
          }
          else
          {
            cout << "nextd: " << nextd << " change lane disabled! current position: " << vx[0] << "," << vy[0] << endl;
          }

          vector<double> localx(next_x_vals.size());
          vector<double> localy(next_x_vals.size());
          for (int i=0; i < next_x_vals.size(); i++)
          {
            float next_x = (next_x_vals[i] - car_x);
            float next_y = (next_y_vals[i] - car_y);
            localx[i] = next_x*cos(angle) + next_y*sin(angle);
            localy[i] = -next_x*sin(angle) + next_y*cos(angle);
          }

          // fit a polynomial
          smooth_local.set_points(localx, localy);

          // calculate the smoother path
          double localxx = 0.;
          double localxy = 0.;
          for(int i = 0; i < num_points; i++)
          {
            localy[i] = smooth_local(localx[i]);
            double dist = distance(localxx, localxy, localx[i], localy[i]);
            if (dist > smooth_speed(double(i)))
            {
               double heading = atan2(localy[i]-localxy,localx[i]-localxx);
               localx[i] = localxx + smooth_speed(double(i))*cos(heading);
               localy[i] = smooth_local(localx[i]);
               dist = distance(localxx, localxy, localx[i], localy[i]);
            }
            localxx = localx[i];
            localxy = localy[i];
          }

          // convert back to global coordinates
          for (int i=0; i<num_points; i++)
          {
            next_x_vals[i] = localx[i]*cos(angle) - localy[i]*sin(angle) + car_x;
            next_y_vals[i] = localx[i]*sin(angle) + localy[i]*cos(angle) + car_y;
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          for (int i=path_size; i<num_points; i++)
          {
            if (i > 0)
              xyd.push_back(distance(next_x_vals[i-1], next_y_vals[i-1], next_x_vals[i], next_y_vals[i]));
            else
              xyd.push_back(0.);
          }

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          this_thread::sleep_for(chrono::milliseconds(500));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          timestep++;        
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

