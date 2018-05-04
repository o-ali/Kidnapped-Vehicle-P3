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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 100;
	default_random_engine gen;
	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);

	for(int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);

		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;


	for(int i = 0; i < num_particles; i++)
	{
		double new_x, new_y, new_theta;

		if(fabs(yaw_rate) == 0)
		{
			new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			new_theta = particles[i].theta;
		}
		else
		{
			new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+(yaw_rate*delta_t))-sin(particles[i].theta));
			new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta + (yaw_rate*delta_t)));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}

		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);

		//adding some noise to particle filter
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
	double distance, landmark;

	for(int i = 0; i < observations.size(); i++)
	{
		double min_d = 9999999;
		int landmark = -1;
		for(int j = 0; j < predicted.size(); j++)
		{
			/*check distance for observations and change landmark to closest*/
			distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if(distance < min_d)
			{
				min_d = distance;
				landmark = predicted[j].id;
			}
		}
		observations[i].id = landmark;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	/* calculate normalization term
	gauss_norm= (1/(2 * np.pi * sig_x * sig_y))
	// calculate exponent
	exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + ((y_obs - mu_y)**2)/(2 * sig_y**2)
	// calculate weight using normalization terms and exponent
	weight= gauss_norm * math.exp(-exponent)*/
	weights.clear();

	for(int p = 0; p < num_particles; p++)
	{

		vector<LandmarkObs> trans_observations;
		LandmarkObs obs;

		//Coordinate transformation vehicle->map
		for(int i = 0; i < observations.size(); i++)
		{
			LandmarkObs trans_obs;
			obs = observations[i];
			trans_obs.x = cos(particles[p].theta)*obs.x - sin(particles[p].theta)*obs.y + particles[p].x;
			trans_obs.y = sin(particles[p].theta)*obs.x + cos(particles[p].theta)*obs.y + particles[p].y;
			trans_obs.id = obs.id;
			trans_observations.push_back(trans_obs);
		}

		//create list of landmarks to predict(within sensor range)
		vector<LandmarkObs> pred_landmarks;

		for(int i = 0; i < map_landmarks.landmark_list.size(); i++)
		{
			double landmark_x = map_landmarks.landmark_list[i].x_f;
			double landmark_y = map_landmarks.landmark_list[i].y_f;
			double landmark_id = map_landmarks.landmark_list[i].id_i;

			if((landmark_x >= particles[p].x-sensor_range) && (landmark_y >= particles[p].y - sensor_range))
			{
				pred_landmarks.push_back(LandmarkObs{landmark_id,landmark_x, landmark_y});
			}
		}

		dataAssociation(pred_landmarks, trans_observations);

		//reset weight to 1
		particles[p].weight = 1.0;
		double prob = 1;
		double pred_x;
		double pred_y;

		//loop on the coordinates and calculate weight for particles
		//that match the landmarks
		for(int i = 0; i < trans_observations.size(); i++)
		{
			double trans_x = trans_observations[i].x;
			double trans_y = trans_observations[i].y;

			for(int j = 0; j < pred_landmarks.size(); j++)
			{
				if(trans_observations[i].id == pred_landmarks[j].id)
				{
					pred_x = pred_landmarks[j].x;
					pred_y = pred_landmarks[j].y;

					double lm_x = std_landmark[0];
					double lm_y = std_landmark[1];
					//multi-variate gaussian
					prob = (1/(2*M_PI*lm_x*lm_y)) * exp(-1 * (pow(trans_x-pred_x,2)/(2*pow(lm_x, 2)) + (pow(trans_y-pred_y,2)/(2*pow(lm_y, 2)))));
					particles[p].weight *= prob;
				}
			}

		}

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<double> weights;
	for(int i = 0; i < num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}

	default_random_engine gen;
	double max_w = *max_element(weights.begin(), weights.end());

	uniform_int_distribution<int> random_index(0,num_particles-1);
	int index = random_index(gen);

	double beta = 0;
	
	vector<Particle> resample_particles;
	uniform_real_distribution<double> random_weight(0,max_w);

	for(int i = 0; i < num_particles; i++)
	{
		beta += 2*random_weight(gen);
		while(beta > weights[index])
		{
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resample_particles.push_back(particles[index]);
	}

	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
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
