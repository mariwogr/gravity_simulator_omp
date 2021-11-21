//
// Created by mariwogr on 16/10/21.
//

#ifndef GRAVITY_SIMULATOR_SIM_AOS_HPP
#define GRAVITY_SIMULATOR_SIM_AOS_HPP

#define INVALID_NUM_ARGS 1
#define INVALID_NUM_OBJS 2
#define INVALID_NUM_ITER 4
#define INVALID_SEED     8
#define INVALID_SIZE     16
#define INVALID_TIME_STP 32
#define G (6.674 * 1E-11)
#define NUM_THREADS 8

struct parameters {
    const int num_objects;
    int num_iterations;
    int random_seed;
    double size_enclosure;
    double time_step;
};

struct set {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double m;
    bool active;
};

int check_bounce(set *objects, int obj, double size);
int check_collision(set *objects, int i, int j);
int collision_objects(set *objects,int i, int j);
int gravitational_force(int num_objects, set *objects, double time_step, double *force, double *accel);
int parser(int argc, char* argv[]);
int print_error_args(int argc, char* argv[], int error_code);
int write_config(int id, parameters system_data, set *params);

#endif //GRAVITY_SIMULATOR_SIM_AOS_HPP
