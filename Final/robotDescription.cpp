#include "robotDescription.h"


extern ofstream logging_output_file;  // defined in main.cpp
/*
 * Dynamic Cantilever Example
 *
 * Define your soft robot structure(s), boundary conditions,
 * custom external forces, and loggers in the function below.
 */

void get_robot_description(int argc, char** argv,
                           const shared_ptr<softRobots>& soft_robots,
                           const shared_ptr<forceContainer>& forces,
                           shared_ptr<worldLogger>& logger,
                           simParams& sim_params) {

    sim_params.dt = 1e-4;
    sim_params.sim_time = 0.1;
    sim_params.dtol = 1e-3;
    sim_params.show_mat_frames = false;
    //sim_params.enable_2d_sim = false;
    sim_params.nis = BACKWARD_EULER;

    int n = 11;
    double radius = 0.003;
    double radius2 = 2*radius;
    double young_mod = 5e9;
    double young_mod2 = 1e6;
    double density = 1200; //in kg/m^3
    double poisson = 0.1;

    // Create a beam along the x-y plane
    double s_factor = 0.3;
    double bowh = 1.0;
    double cross = 0.05;
    double ccross = 2*cross;
    double midbow = bowh/2;
    double band = 0.15;
    soft_robots->addLimb(Vector3d(0.0, 0.0, 0.0)*s_factor, Vector3d(ccross, 0.0, 0.0)*s_factor, n, density, radius, young_mod, poisson); //limb0 bottom cross
    soft_robots->addLimb(Vector3d(0.0, 0.0, 0.0)*s_factor, Vector3d(0.0, 0.0, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb1 first bow
    soft_robots->addLimb(Vector3d(cross, -cross, 0.0)*s_factor, Vector3d(cross, 0.0, 0.0)*s_factor, n, density, radius, young_mod, poisson); //limb2 bottom middle cross
    soft_robots->addLimb(Vector3d(cross, cross, 0.0)*s_factor, Vector3d(cross, 0.0, 0.0)*s_factor, n, density, radius, young_mod, poisson); //limb3 bottom middle cross
    soft_robots->addLimb(Vector3d(ccross, 0.0, 0.0)*s_factor, Vector3d(ccross, 0.0, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb4 second bow
    soft_robots->addLimb(Vector3d(cross, -cross, 0.0)*s_factor, Vector3d(cross, -cross, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb5 third bow
    soft_robots->addLimb(Vector3d(cross, cross, 0.0)*s_factor, Vector3d(cross, cross, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb6 fourth bow
    soft_robots->addLimb(Vector3d(0.0, 0.0, bowh)*s_factor, Vector3d(ccross, 0.0, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb7 top center cross
    soft_robots->addLimb(Vector3d(cross, -cross, bowh)*s_factor, Vector3d(cross, 0.0, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb8 top cross
    soft_robots->addLimb(Vector3d(cross, cross, bowh)*s_factor, Vector3d(cross, 0.0, bowh)*s_factor, n, density, radius, young_mod, poisson); //limb9 top cross
    /*
    soft_robots->addLimb(Vector3d(cross, 0.0, 0.0)*s_factor, Vector3d(cross, 0.0, midbow)*s_factor, n, density, radius, young_mod2/100, poisson); //limb10 center rod
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(0.0, 0.0, midbow+band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb11 elastic band1 first bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(0.0, 0.0, midbow-band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb12 elastic band2 first bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(ccross, 0.0, midbow+band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb13 elastic band3 second bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(ccross, 0.0, midbow-band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb14 elastic band4 second bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(cross, -cross, midbow+band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb15 elastic band5 third bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(cross, -cross, midbow-band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb16 elastic band6 third bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(cross, cross, midbow+band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb17 elastic band7 fourth bow
    soft_robots->addLimb(Vector3d(cross, 0.0, midbow)*s_factor, Vector3d(cross, cross, midbow-band)*s_factor, n, density, radius2, young_mod2, 0.5); //limb18 elastic band8 fourth bow
*/
    
    //joint on which limb, joint on which node
    soft_robots->createJoint(0,0); // joint 0 bottom cross to first bow
    soft_robots->createJoint(0,round(n/2)); // joint 1 bottom cross middle
    soft_robots->createJoint(0,-1); // joint 2 bottom cross to second bow
    soft_robots->createJoint(1,-1); // joint 3 first bow to top cross
    soft_robots->createJoint(2,0); // joint 4 bottom cross to third bow
    soft_robots->createJoint(3,0); // joint 5 bottom cross to fourth bow
    soft_robots->createJoint(4,-1); // joint 6 second bow to top cross
    soft_robots->createJoint(5,-1); // joint 7 third bow to top cross
    soft_robots->createJoint(6,-1); // joint 8 fourth bow to top cross
    soft_robots->createJoint(7,round(n/2)); // joint 9 top cross middle
    /*
    soft_robots->createJoint(10,-1); // joint 10 center rod top
    soft_robots->createJoint(1,round(3*n/4)); // joint 11 elastic band1 to first bow
    soft_robots->createJoint(1,round(n/4)); // joint 12 elastic band2 to first bow
    soft_robots->createJoint(4,round(3*n/4)); // joint 13 elastic band3 to second bow
    soft_robots->createJoint(4,round(n/4)); // joint 14 elastic band4 to second bow
    soft_robots->createJoint(5,round(3*n/4)); // joint 15 elastic band5 to third bow
    soft_robots->createJoint(5,round(n/4)); // joint 16 elastic band6 to third bow
    soft_robots->createJoint(6,round(3*n/4)); // joint 17 elastic band7 to fourth bow
    soft_robots->createJoint(6,round(n/4)); // joint 18 elastic band8 to fourth bow
    */

    //joint #, limb connect to the joint, node of the limb
    soft_robots->addToJoint(0,1,0); // bottom cross to first bow
    soft_robots->addToJoint(1,2,-1); // bottom cross middle
    soft_robots->addToJoint(1,3,-1); // bottom cross middle
    soft_robots->addToJoint(2,4,0); // bottom cross to second bow
    soft_robots->addToJoint(4,5,0); // bottom cross to third bow
    soft_robots->addToJoint(5,6,0); // bottom cross to fourth bow
    soft_robots->addToJoint(3,7,0); // first bow to top cross
    soft_robots->addToJoint(6,7,-1); // second bow to top cross
    soft_robots->addToJoint(7,8,0); // third bow to top cross
    soft_robots->addToJoint(8,9,0); // fourth bow to top cross
    soft_robots->addToJoint(9,8,-1); // top cross middle
    soft_robots->addToJoint(9,9,-1); // top cross middle
    /*
    soft_robots->addToJoint(1,10,0); // bottom cross middle to center rod
    soft_robots->addToJoint(10,11,0); // center rod to elastic band1
    soft_robots->addToJoint(11,11,-1); // elastic band1 to first bow
    soft_robots->addToJoint(10,12,0); // center rod to elastic band2
    soft_robots->addToJoint(12,12,-1); // elastic band2 to first bow
    soft_robots->addToJoint(10,13,0); // center rod to elastic band3
    soft_robots->addToJoint(13,13,-1); // elastic band3 to second bow
    soft_robots->addToJoint(10,14,0); // center rod to elastic band4
    soft_robots->addToJoint(14,14,-1); // elastic band4 to second bow
    soft_robots->addToJoint(10,15,0); // center rod to elastic band5
    soft_robots->addToJoint(15,15,-1); // elastic band5 to third bow
    soft_robots->addToJoint(10,16,0); // center rod to elastic band6
    soft_robots->addToJoint(16,16,-1); // elastic band6 to third bow
    soft_robots->addToJoint(10,17,0); // center rod to elastic band7
    soft_robots->addToJoint(17,17,-1); // elastic band7 to fourth bow
    soft_robots->addToJoint(10,18,0); // center rod to elastic band8
    soft_robots->addToJoint(18,18,-1); // elastic band8 to fourth bow
    */
    
    soft_robots->setup(); 
    // to run the codes:
    // cp examples/cantilever_case/cantileverExample.cpp robotDescription.cpp [this read the cantilever example and overwrite! the values in the example to robotDescription.cpp]
    // [omit the above step and just edit robotDescription for our own project!]

    // cd build && make -j4 && cd .. [build dismech.sh, current directory must be in the main directory]
    // ./dismech.sh [run]

    // Add weight to system
    Vector3d gravity_vec(0.0, 0.0, -9.8);
    forces->addForce(make_shared<gravityForce>(soft_robots, gravity_vec));

    // Create an external constant uniform force
    shared_ptr<uniformConstantForce> uniform_force = make_shared<uniformConstantForce>(soft_robots);
    Vector3d force = Vector3d::Zero();
    string force_value = "165";
    force(2) = stod(force_value);  // apply along z-direction
    // Add uniform constant force to the cylinder
    //uniform_force->add_force_to_limb(0, 2*force);
    //uniform_force->add_force_to_limb(2, 1*force);
    //uniform_force->add_force_to_limb(3, 1*force);

    uniform_force->add_force_to_limb(7, 2*force);
    uniform_force->add_force_to_limb(8, 1*force);
    uniform_force->add_force_to_limb(9, 1*force);
    forces->addForce(uniform_force);

    // Add viscous damping
    //double viscosity = 7.0;
    //forces->addForce(make_shared<dampingForce>(soft_robots, viscosity));

    // Add ground contact

    double delta = 5e-4;
    double nu = 5e-3;
    double mu = 0.0;
    double floor_z = -0.003;
    forces->addForce(make_shared<floorContactForce>(soft_robots, delta, nu, mu, floor_z));


    // Fix one end
    /*
    soft_robots->lockEdge(0,round(n/2));
    soft_robots->lockEdge(0,0);
    soft_robots->lockEdge(0,n);
    soft_robots->lockEdge(2,0);
    soft_robots->lockEdge(3,0);
    */

    // Read initial velocity values
    //vector<Vector3d> velocities;
    //load_txt("examples/cantilever_case/cantilever_init_velocity_n=201.txt", velocities);

    // Apply the velocities
    //soft_robots->applyInitialVelocities(0, velocities);

    // Set logger to record nodes
    string logfile_base = "log_files/project";
    int logging_period = 50;
    logger = make_shared<rodNodeLogger>(logfile_base, logging_output_file, logging_period);

    
}
