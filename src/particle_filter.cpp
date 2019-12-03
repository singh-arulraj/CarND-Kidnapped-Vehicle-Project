/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
#define NUM_PARTICLES 1000
double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  int particle_index;
  std::default_random_engine gen;
  num_particles = NUM_PARTICLES;  // TODO: Set the number of particles
   // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (particle_index = 0; particle_index < num_particles; particle_index++) {
	Particle particle = {.id = particle_index,
                             .x = dist_x(gen),
		             .y = dist_y(gen),
                             .theta = dist_theta(gen),
                             .weight = 1
                            };
        weights.push_back(1);
        particles.push_back(particle);
  }
  is_initialized = true;  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   int particle_id;
   std::default_random_engine gen;
   normal_distribution<double> coord_y(0.0, std_pos[1]);
   normal_distribution<double> coord_x(0.0, std_pos[0]);
   normal_distribution<double> coord_theta(0.0, std_pos[2]);
   
   if (fabs(yaw_rate) < 0.0001) {
       yaw_rate = 0.0001;
   }
   

   for (particle_id = 0; particle_id < num_particles; particle_id++) {
        particles[particle_id].x     += (velocity * (sin(particles[particle_id].theta + yaw_rate * delta_t) - sin(particles[particle_id].theta)) / yaw_rate);
        particles[particle_id].y     += (velocity * (cos(particles[particle_id].theta) - cos(particles[particle_id].theta + yaw_rate * delta_t)) / yaw_rate);
	particles[particle_id].theta += (delta_t * yaw_rate);
        //Adding Gaussian Noise
        /*particles[particle_id].x +=  coord_x(gen);
        particles[particle_id].y +=  coord_y(gen);
        particles[particle_id].theta += coord_theta(gen); 
   */
}
   
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   int numObservations = observations.size();
   int numParticles = particles.size();
   float distance, temp_dist, weight;
   int index;
   float sum_weights = 0;
   for (int particleIndex = 0; particleIndex < numParticles; particleIndex++) {
   	Particle &particle = particles[particleIndex];
   	float obs_x ;
   	float obs_y;
   	particle.weight = 1;
        for (int i= 0; i < observations.size(); i++) {
   		//Convert Observation point to map Coordinates
		obs_x = observations[i].x * cos(particle.theta) - observations[i].y * sin(particle.theta) + particle.x;
   		obs_y = observations[i].x * sin(particle.theta) + observations[i].y * cos(particle.theta) + particle.y;
                
		//Generate the associations
                distance = -1;       
                for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			temp_dist = dist (obs_x,obs_y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
			if (temp_dist > sensor_range) {
				continue;
			}
			if (distance == -1) {
                            distance = temp_dist;
                            index = j;
                        } else if (distance > temp_dist) {
				distance = temp_dist;
				index = j;
			}
		}
                if (distance != -1) {
		weight = multiv_prob(std_landmark[0], std_landmark[1], particle.x, particle.y, map_landmarks.landmark_list[index].x_f, map_landmarks.landmark_list[index].y_f);
		particle.weight *= weight;
                } else {
particle.weight = 0;
}
   	}
        sum_weights += particle.weight;
        
   }
  /*
 for (int particleIndex = 0; particleIndex < numParticles; particleIndex++) {
       particles[particleIndex].weight /= sum_weights;
}
*/	

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   	// Update weights`
	std::default_random_engine gen;
        weights.clear();
        int numParticles = particles.size();
	for(int i=0; i<numParticles; ++i){
		weights.push_back(particles[i].weight);
	}

	std::discrete_distribution<int> particle_dist(weights.begin(),weights.end());

	// Resample particles
	vector<Particle> new_particles;
	new_particles.resize(num_particles);
	for(int i=0; i<num_particles; ++i){

		auto index = particle_dist(gen);
		new_particles[i] = particles[index];
	}
	particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
