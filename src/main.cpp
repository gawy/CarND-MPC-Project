#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

// for convenience
using json = nlohmann::json;

const double LATENCY = 0.1; // 100 ms
const double Lf = 2.67;

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
  auto b2 = s.rfind("}]");
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

void convertMapToCarCoords(vector<double> &ptsx, vector<double> &ptsy,
                      vector<double> &outx, vector<double> &outy,
                      double xCenter, double yCenter, double psi) {
  for (int i = 0; i < ptsx.size(); ++i) {
    double x0 = ptsx[i];
    double y0 = ptsy[i];
    double x = -(x0 - xCenter) * cos(M_PI - psi) + (y0 - yCenter) * sin(M_PI - psi);
    double y = -(x0 - xCenter) * sin(M_PI - psi) - (y0 - yCenter) * cos(M_PI - psi);

    outx[i] = x;
    outy[i] = y;
  }
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

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          milliseconds time_start = duration_cast< milliseconds >(system_clock::now().time_since_epoch());


          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double acceleration = j[1]["throttle"];
          double delta = -1.0 * (double)j[1]["steering_angle"];

          cout << "Steering angle: " << j[1]["steering_angle"] << ", converted=" << delta << endl;

          // convert speed to m/s
          v = v * 1609.34 / 3600;

          /*
          * Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */

          double x_car = px;
          double y_car = py;
          double psi_car = psi;

          vector<double> car_x(ptsx.size());
          vector<double> car_y(ptsx.size());

          convertMapToCarCoords(ptsx, ptsy, car_x, car_y, x_car, y_car, psi_car);
          Eigen::Map<VectorXd> x_pts(car_x.data(), car_x.size());
          Eigen::Map<VectorXd> y_pts(car_y.data(), car_y.size());


          x_car = 0.0;
          y_car = 0.0;
          psi_car = 0.0;

          cout << "X size: "<< x_pts.size()<< endl;
          VectorXd coeffs = polyfit(x_pts, y_pts, 3); //fit polinomial to waypoints


          cout << "Poly test: x=" << x_pts[1] << ", y=" << polyeval(coeffs, x_pts[1])
               << "; x1=" << x_pts[4] << ", y1=" << polyeval(coeffs, x_pts[4])
               << "; x0=" << 0 << ", y0=" << polyeval(coeffs, 0)
               << endl;

          double cte = y_car - polyeval(coeffs, x_car);
          double tan_psi = coeffs[1] + 2*coeffs[2]*x_car + 3*coeffs[3]*x_car*x_car;

          double epsi = psi_car - atan(tan_psi);

          cout << "epsi=" << epsi << ", k="<<tan_psi << "(" << rad2deg(tan_psi) << " deg) " << endl;

          if (LATENCY > 0.0) {
            //adjust for latency
            // simple way is to move car position to after the latency time
            x_car += v * cos(psi_car) * LATENCY;
            y_car += v * sin(psi_car) * LATENCY;
            psi_car += v / Lf * delta * LATENCY;

            cout << "Initial state corrected for Latency" << endl;
            cout << "x=" << x_car << ", y=" << y_car << ", psi=" << psi_car << ", v=" << v << endl;

            cte = cte + v * sin(epsi) * LATENCY;
            epsi = epsi + v / Lf * delta * LATENCY;

            v += acceleration * LATENCY;

            px = x_car; py = y_car; psi = psi_car;
          }

          Eigen::VectorXd state(6);
          state << x_car, y_car, psi_car, v, cte, epsi;

          cout << "Polynomial coeffs: " << coeffs[0] << "; " << coeffs[1]<< "; " << coeffs[2]<< "; " << coeffs[3] << endl;

          long pts_len = x_pts.size();
          vector<double> controls = mpc.Solve(state, coeffs, x_pts[pts_len - 1], y_pts[pts_len - 1]);

          double steer_value = controls[0] / deg2rad(25);
          double throttle_value = controls[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          for (int i = 0; i < mpc.x_waypts.size(); ++i) {
            mpc_x_vals.push_back(mpc.x_waypts[i]);
            mpc_y_vals.push_back(mpc.y_waypts[i]);
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals(ptsx.size());
          vector<double> next_y_vals(ptsy.size());

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          for (int k = 0; k < x_pts.size(); ++k) {
            next_x_vals[k] = x_pts[k];
            next_y_vals[k] = y_pts[k];
          }


          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          milliseconds time_end = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
          auto processing_time = (time_end - time_start).count();
          cout << "Processing took " << processing_time << "ms" << endl;
          if (processing_time > 100) {
            cout << "!!! Took more than 100ms for MPC solver" << endl;
          }

          if (LATENCY > 0.0) {
            this_thread::sleep_for(chrono::milliseconds((int) (LATENCY * 1000)));
          }

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

