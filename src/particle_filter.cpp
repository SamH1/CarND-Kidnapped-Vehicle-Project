/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

//output logfile for debugging
string out_file_name_ = "../LogFile.txt";
//std::ofstream out_file_(out_file_name_.c_str(), std::ofstream::out);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	//set number of particles
	num_particles = 100;

	//set random generator engine
	default_random_engine gen;

	//set initial Gaussian noise
	normal_distribution<double> N_x(x, std[0]);           
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);

	//initialize particles
	for (int i = 0; i < num_particles; i++){
		Particle p;
		p.id = i;
		p.x = N_x(gen);
		p.y = N_y(gen);
		p.theta = N_theta(gen);
		p.weight = 1;
		particles.push_back(p);
		weights.push_back(1);
	}

	//particle filter initialized
	is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for (int i = 0; i < num_particles; i++){
		double x1;
		double y1;
		double theta1;

		//updating x, y and the yaw : Lesson 14 - section 6 to 8

		Particle p_i = particles[i];

		// if yaw_rate ==0 calculate traveled distance and project on x and y axis...
		if (yaw_rate == 0){
			x1 = p_i.x + velocity*delta_t*cos(p_i.theta);
			y1 = p_i.y + velocity*delta_t*sin(p_i.theta);
			theta1 = p_i.theta;

		} else {
			// include yaw movement (Lesson 14 - section 8)
			x1 = p_i.x + (velocity/yaw_rate) * (sin(p_i.theta + delta_t*yaw_rate) - sin(p_i.theta));
			y1 = p_i.y + (velocity/yaw_rate) * (cos(p_i.theta) - cos(p_i.theta + delta_t*yaw_rate));
			theta1 = p_i.theta + delta_t*yaw_rate;			
		}

		// add Gaussian noise to predicted positions
		normal_distribution<double> N_x(x1, std_pos[0]);           
		normal_distribution<double> N_y(y1, std_pos[1]);
		normal_distribution<double> N_theta(theta1, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// loop over all observations and landmarks and store the landmak ID with min distance	
	for (int i = 0; i < observations.size(); i++){
		double min_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
		int landmark_id = 0;

		for (int j= 0; j < predicted.size(); j++){

			double dist_ = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

			if (dist_  < min_dist){
				landmark_id = j;
				min_dist = dist_;
			}
		}

		observations[i].id = landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	///// TBC not sure if I should start with this loop
	for (int p = 0; p < particles.size(); p++){
		/// peut être qu'il manque qq chose ici aussi une initialisation ou qq chose...

		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		vector<LandmarkObs> trans_observations;
		LandmarkObs obs;
		
		// Transform observation from car to map coordinates
		for (int i = 0; i < observations.size(); i++){
			LandmarkObs trans_obs;
			obs = observations[i];

			//transformation from vehicule to map coordinates: rotation based on the particle theta and then x, y translation
			trans_obs.x = particles[p].x + obs.x*cos(particles[p].theta) - obs.y*sin(particles[p].theta);
			trans_obs.y = particles[p].y + obs.x*sin(particles[p].theta) + obs.y*cos(particles[p].theta);
			trans_observations.push_back(trans_obs);
		}

		// Find landmarks around the particle within LiDAR range
		vector<LandmarkObs> landmarks_inRange;		
		for (int j = 0;  j < map_landmarks.landmark_list.size(); j++) {
			int id_ = map_landmarks.landmark_list[j].id_i;
			double x_ = map_landmarks.landmark_list[j].x_f;
			double y_ = map_landmarks.landmark_list[j].y_f;

			double dist_ =  dist(x_, y_, particles[p].x, particles[p].y);

			if (dist_ < sensor_range) {
				LandmarkObs landmark_ = {id_, x_, y_};
				landmarks_inRange.push_back(landmark_);
			}
		}

		//associate obs with landmarks ids
		dataAssociation(landmarks_inRange, trans_observations);

		//initialise particle weight
		particles[p].weight = 1.0;
		
		//update weight
		for (int i = 0; i < trans_observations.size(); i++){
			if (trans_observations[i].id != 0){
				double meas_x = trans_observations[i].x;
				double meas_y = trans_observations[i].y;
				double mu_x = landmarks_inRange[trans_observations[i].id].x;
				double mu_y = landmarks_inRange[trans_observations[i].id].y;
	/*
				out_file_ << "x1 " << meas_x << endl;
				out_file_ << "y1 " << meas_y << endl;
				out_file_ << "x2 " << mu_x << endl;
				out_file_ << "y2 " << mu_y << endl;
				out_file_ << "std_x " << std_landmark[0] << endl;
				out_file_ << "std_y " << std_landmark[1] << endl;
	*/
				// Section 13 - Calculating the Particle's Final Weight
				double A = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);
				double B = pow(meas_x-mu_x, 2)/(2*pow(std_landmark[0],2));
				double C = pow(meas_y-mu_y, 2)/(2*pow(std_landmark[1],2));

				double P = A * exp(-(B+C));
			
			
	//			out_file_ << "P " << P << endl;
	//			out_file_ << "-------------------------------------------" << endl;
				if(P > 0){
					particles[p].weight *= P;
				}
				
				associations.push_back(trans_observations[i].id);
				sense_x.push_back(trans_observations[i].x);
				sense_y.push_back(trans_observations[i].y);
			}
		}
		particles[p] = SetAssociations(particles[p], associations, sense_x, sense_y);
		weights[p] = particles[p].weight;
	}
}  

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resampled_particles;

	for (int i = 0; i < num_particles; i++){
		resampled_particles.push_back(particles[distribution(gen)]);
	}

	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
