#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <omp.h>

#include "psoa.hpp"

using namespace std;

/* *
 * This function will check the parameters at the beginning
 *
 * @param int argc                 its the number about how many parameters are in the execution
 * @param char argv             it's an array of chars, inside it, we have the arguments
 * @return 0 on success
 */

int parser(int argc, char* argv[]){
    int ret = 0;
    //Checking if the number of arguments its correct
    if (argc != 6 ){
        ret |= INVALID_NUM_ARGS;
    } else{
	    //checking if the number of objs is smaller than zero
	    if ( stoi(argv[1]) <= 0 )
		    ret |= INVALID_NUM_OBJS;

	    //checking if the number of iterations is smaller than zero
	    if ( stoi(argv[2]) < 0 )
		    ret |= INVALID_NUM_ITER;

	    //checking the seed if it's a positive number
	    if ( stoi(argv[3]) <= 0 )
		    ret |= INVALID_SEED;

	    //checking if size_enclosure is positive
	    if ( stod(argv[4]) <= 0.0 )
		    ret |= INVALID_SIZE;

	    //checking if time_step is a real number positive
	    if ( stod(argv[5]) <= 0.0 )
		    ret |= INVALID_TIME_STP;
    }
    return ret;
}


/*
 * This function updates the speed vector v and the position of every point in the set objects of points
 *
 * @param: int num_objects          the total of points in the set of objects
 * @param: set objects              structure of objects, with all the components of each point in the simulator
 * @param: float time_step          time step to obtain the speed and position of the point
 *
 * @return 0                        if the function was executed correctly
 */

int gravitational_force(int num_objects, set objects, double time_step, double *force, double *accel) {

    double powSqX;
    double powSqY;
    double powSqZ;
    double norm;
    double fx;
    double fy;
    double fz;
    
    // The execution will pass through two nested loops to obtain the sum of gravitational forces of every point with
    // the other points. Analogous to take a screenshot of the system before updating speeds and positions.
    #pragma omp parallel
    {
        double powSqX;
        double powSqY;
        double powSqZ;
        double norm;
        double fx;
        double fy;
        double fz;

        int t = omp_get_thread_num();

        for (int i = t; i < num_objects; i = i + NUM_THREADS) {
            if (objects.active[i]) {
                for (int j = i + 1; j < num_objects; j++) {
                    if (objects.active[j]) {
                        powSqX = (objects.x[j] - objects.x[i]) * (objects.x[j] - objects.x[i]);
                        powSqY = (objects.y[j] - objects.y[i]) * (objects.y[j] - objects.y[i]);
                        powSqZ = (objects.z[j] - objects.z[i]) * (objects.z[j] - objects.z[i]);
                        norm = std::sqrt(powSqX + powSqY + powSqZ);
                        // It will return the three components of the gravitational force between i and j

                        fx = (G * objects.m[i] * objects.m[j] * (objects.x[j] - objects.x[i])) / (norm * norm * norm);
                        fy = (G * objects.m[i] * objects.m[j] * (objects.y[j] - objects.y[i])) / (norm * norm * norm);
                        fz = (G * objects.m[i] * objects.m[j] * (objects.z[j] - objects.z[i])) / (norm * norm * norm);

                        force[3 * i] += fx;
                        force[3 * i + 1] += fy;
                        force[3 * i + 2] += fz;
                        force[3 * j] -= fx;
                        force[3 * j + 1] -= fy;
                        force[3 * j + 2] -= fz;
                    }
                }
            }
        }
    }
    // Once we have a screenshot of the system in force array, update each active object
    for (int i = 0; i < num_objects; i++) {
        if(objects.active[i]) {
            // Updates the acceleration
            accel[0] = 1.0/objects.m[i] * force[i * 3];
            accel[1] = 1.0/objects.m[i] * force[(i * 3) + 1];
            accel[2] = 1.0/objects.m[i] * force[(i * 3) + 2];

            // Updates the speed
            objects.vx[i] = objects.vx[i] + accel[0] * time_step;
            objects.vy[i] = objects.vy[i] + accel[1] * time_step;
            objects.vz[i] = objects.vz[i] + accel[2] * time_step;

            // Updates the position
            objects.x[i] = objects.x[i] + objects.vx[i] * time_step;
            objects.y[i] = objects.y[i] + objects.vy[i] * time_step;
            objects.z[i] = objects.z[i] + objects.vz[i] * time_step;
        }
    }
    return 0;
}

/*
 * This function check if the object bounce with a wall and change the values if it's necessary
 *
 * @param: set objects              structure of objects, with all the components of each point in the simulator
 * @param: float size         It's the size of the wall
 * @param: int obj          it's the object which we are going to check
 *
 * @return 0                        if the function was executed correctly
 */
int check_bounce(set objects, int obj, double size){

    //check if the object bounce with a wall

    if(objects.x[obj] <= 0){
        objects.x[obj] = 0;
        objects.vx[obj] = -1 * objects.vx[obj];
    }

    if(objects.y[obj] <= 0){
        objects.y[obj] = 0;
        objects.vy[obj] = -1 * objects.vy[obj];
    }

    if(objects.z[obj] <= 0){
        objects.z[obj] = 0;
        objects.vz[obj] = -1 * objects.vz[obj];
    }

    if(objects.x[obj] >= size){
        objects.x[obj] = size;
        objects.vx[obj] = -1 * objects.vx[obj];
    }

    if(objects.y[obj] >= size){
        objects.y[obj] = size;
        objects.vy[obj] = -1 * objects.vy[obj];
    }

    if(objects.z[obj] >= size){
        objects.z[obj] = size;
        objects.vz[obj] = -1 * objects.vz[obj];
    }

    return 0;
}

/*
* This function will check if the object collisions with another
*
* @param: set objects          array of objects with their properties
* @param: int i                array position of the first object
* @param: int j                array position of the second object
*/

int check_collision(set objects, int i, int j){
    double distance = std::sqrt((objects.x[i] - objects.x[j]) * (objects.x[i] - objects.x[j]) \
                            + (objects.y[i] -objects.y[j]) * (objects.y[i] -objects.y[j]) \
                            + (objects.z[i] -objects.z[j]) * (objects.z[i] -objects.z[j]));

    if(distance < 1.0){
        collision_objects(objects, i, j);
    }
    return 0;
}


/*
* This function will update the objects and their collisions
*
* @param: set objects          array of objects with their properties
* @param: int i                array position of the first object
* @param: int j                array position of the second object
*/
int collision_objects(set objects, int i, int j){
    objects.m[i] += objects.m[j];
    objects.vx[i] += objects.vx[j];
    objects.vy[i] += objects.vy[j];
    objects.vz[i] += objects.vz[j];

    objects.active[j] = false;
    return 0;
}

/* *
 * This function will write the errors in the parameter in error case
 *
 * @param int argc                 This is the number of arguments
 * @param char* argv is a pointer to array of chars (strings) which it has the values
 *
 * @return 0 on success
 */
int print_error_args(int argc, char* argv[], int error_code) {
    /*This function will print in the standard output the parameters when the function was called
      and it will show the errors while doing it.*/
    if(error_code & INVALID_NUM_OBJS)
    	cerr << "Error: Invalid number of objects" << endl;
    if(error_code & INVALID_NUM_ITER)
    	cerr << "Error: Invalid number of iterations" << endl;
    if(error_code & INVALID_SEED)
    	cerr << "Error: Invalid seed" << endl;
    if(error_code & INVALID_SIZE)
    	cerr << "Error: Invalid size enclosure" << endl;
    if(error_code & INVALID_TIME_STP)
    	cerr << "Error: Invalid time step" << endl;
    if(error_code & INVALID_NUM_ARGS)
    	cerr << "Error: Wrong number of parameters" << endl;
    cerr << argv[0] << " invoked with " << argc - 1  << " parameters." << endl;
    cerr << "Arguments:" << endl;

    /* If argc is not 5, it will show a ? character when a parameter is not in the function call */
    if (1 < argc) { cerr << "  num_objects: " << argv[1] << endl; }
    else { cerr << "  num_objects: ?" << endl; }

    if (2 < argc) { cerr << "  num_iterations: " << argv[2] << endl; }
    else { cerr << "  num_iterations: ?" << endl; }

    if (3 < argc) { cerr << "  random_seed: " << argv[3] << endl; }
    else { cerr << "  random_seed: ?" << endl; }

    if (4 < argc) { cerr << "  size_enclosure: " << argv[4] << endl; }
    else { cerr << "  size_enclosure: ?" << endl; }

    if (5 < argc) { cerr << "  time_step: " << argv[5] << endl; }
    else { cerr << "  time_step: ?" << endl; }

    return 0;
}

/* *
 * This function will write the parameters from the main program in the init_config file or in the final_config file
 *
 * @param int id                 whether it's the first or last file
 * @param parameters system_data data of the systema (size_enclosure, etc.)
 * @param set objects            structure containing the information of the objects
 * @return 0 on success
 */
int write_config(int id, parameters system_data, set objects){
    ofstream out_file;
    char res[5001];

    /*If the id is 0 it will write the content in the init_config file*/
    if (id == 0){
        out_file.open("init_config.txt");
    }
        /*If the id is different from 0 the content will be written in the final_config file*/
    else { out_file.open("final_config.txt"); }

    sprintf(res, "%.3f ", system_data.size_enclosure);
    out_file << res;
    sprintf(res, "%.3f ", system_data.time_step);
    out_file << res;
    sprintf(res, "%d", system_data.num_objects);
    out_file << res << endl;

    for(int i = 0; i < system_data.num_objects; i++){
        if(objects.active[i]) {
            sprintf(res,
                    "%.3f %.3f %.3f %.3f %.3f %.3f %.3f",
                    objects.x[i], objects.y[i], objects.z[i], objects.vx[i], objects.vy[i], objects.vz[i],
                    objects.m[i]);
            out_file << res << endl;
        }

    }
    out_file.close();
    return 0;
}

int main(int argc, char* argv[]) {
    /*The array of parameters argv passes through a parser to check all the arguments are correct*/
    int retcode = parser(argc, argv);
    /*The result of the parser will be equal to -1 or -2 if there are errors with the arguments*/
    if(retcode != 0){
        /*If there is an error, the program will call the function print_error_args to print them
         * through the standard output*/
        print_error_args(argc, argv, retcode);

        /*The main function will return retcode, which can be equal to -1 if there are not enough
         * arguments and -2 if the arguments are not correct*/
        if (retcode == 1){
        	return -1;
        } else{
        	return -2;
        }
    }

    /* Store simulation arguments in a structure */
    parameters system_data{ stoi(argv[1]), stoi(argv[2]),
                            stoi(argv[3]), stod(argv[4]),
                            stod(argv[5])};

    /* Declare the structure that holds objects' information */
    set objects{
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (double *) malloc(sizeof(double) * system_data.num_objects),
            (bool *) malloc(sizeof(bool) * system_data.num_objects)
    };

    cout << "Creating simulation:" << endl;
    cout << "  num_objects: " << system_data.num_objects << endl;
    cout << "  num_iterations: " << system_data.num_iterations << endl;
    cout << "  random_seed: " << system_data.random_seed << endl;
    cout << "  size_enclosure: " << system_data.size_enclosure << endl;
    cout << "  time_step: " << system_data.time_step << endl;

    /* Create mersenne-twister generator and create a uniform and a normal distribution */
    mt19937_64 gen64(system_data.random_seed);
    uniform_real_distribution<> position_unif_dist(0, system_data.size_enclosure);
    normal_distribution<> mass_norm_dist{1E21, 1E15};

    double *force = (double *) malloc(sizeof(double) * system_data.num_objects * 3);
    double accel[3] = {0,0,0};

    /* Initialize x, y, z and m attributes of each object */
    for(int i = 0; i < system_data.num_objects; i++){
            objects.x[i] = position_unif_dist(gen64);
            objects.y[i] = position_unif_dist(gen64);
            objects.z[i] = position_unif_dist(gen64);
            objects.m[i] = mass_norm_dist(gen64);
            objects.active[i] = true;
    }

    /* Write initial configuration to a file*/
    write_config(0, system_data, objects);


    /* Initial collision checking */
    for(int i = 0; i < system_data.num_objects; i++){
        if( !objects.active[i] ){ continue; }
        for(int j = i + 1; j < system_data.num_objects; j++){
            if(objects.active[j])
                check_collision(objects, i, j);
        }
    }

    /* Body of the simulation */
    for(int i = 0; i < system_data.num_iterations; i++){
        for(int foo=0; foo < system_data.num_objects * 3; foo++){force[foo] = 0;}
            gravitational_force(system_data.num_objects, objects, system_data.time_step, force, accel);

            for(int a = 0; a < system_data.num_objects; a++){
                if (objects.active[a])
                    check_bounce(objects, a, system_data.size_enclosure);
            }
            for(int a = 0; a < system_data.num_objects; a++){
                if(objects.active[a]){
                    for(int b = a + 1; b < system_data.num_objects; b++){
                        if (objects.active[b]){
                            check_collision(objects, a, b);
                    }
                }
            }
        }
    }

    /* Write final configuration to a file */
    write_config(1, system_data, objects);
    free(force);
    return 0;
}
