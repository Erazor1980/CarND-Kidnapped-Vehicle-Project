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

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 50;

	/* store standard deviations for easier usage */
	const double std_x		= std[ 0 ];
	const double std_y		= std[ 1 ];
	const double std_yaw	= std[ 2 ];

	/* creating normal Gaussian distributions */
	default_random_engine gen;
	normal_distribution< double > dist_x( x, std_x );
	normal_distribution< double > dist_y( y, std_y );
	normal_distribution< double > dist_theta( theta, std_yaw );

	for( int i = 0; i < num_particles; ++i )
	{
		Particle newParticle;
		newParticle.x		= dist_x( gen );
		newParticle.y		= dist_y( gen );
		newParticle.theta	= dist_theta( gen );
		newParticle.id		= i;
		newParticle.weight	= 1.0;

		particles.push_back( newParticle );
		weights.push_back( 1.0 );
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	/* apply motion models to all particles */
	for( int i = 0; i < particles.size(); ++i )
	{
		Particle& p = particles[ i ];

		/* two different motion models, depending on yaw_rate! */
		if( abs( yaw_rate ) < 0.000001 )		// yaw rate == 0
		{
			p.x = p.x + velocity * delta_t * cos( yaw_rate );
			p.y = p.y + velocity * delta_t * sin( yaw_rate );
			p.theta = 0.0;
		}
		else                                    // yaw rate != 0
		{
			p.x = p.x + velocity / yaw_rate * ( sin( p.theta + yaw_rate * delta_t ) - sin( p.theta ) );
			p.y = p.y + velocity / yaw_rate * ( cos( p.theta ) - cos( p.theta + yaw_rate * delta_t ) );
			p.theta = p.theta + yaw_rate * delta_t;
		}
	}


	// add random Gaussian noise to all particles
	/* store standard deviations for easier usage */
#if 0 // old version
	const double std_x		= std_pos[ 0 ];
	const double std_y		= std_pos[ 1 ];
	const double std_yaw	= std_pos[ 2 ];

	for( int i = 0; i < particles.size(); ++i )
	{
		Particle& p = particles[ i ];

		/* creating normal Gaussian distributions */
		default_random_engine gen;
		normal_distribution< double > dist_x( p.x, std_x );
		normal_distribution< double > dist_y( p.y, std_y );
		normal_distribution< double > dist_theta( p.theta, std_yaw );

		p.x = dist_x( gen );
		p.y = dist_y( gen );
		p.theta = dist_theta( gen );
	}
#else   // new version 
	const double std_x = std_pos[ 0 ];
	const double std_y = std_pos[ 1 ];
	const double std_yaw = std_pos[ 2 ];

	default_random_engine gen;
	normal_distribution< double > dist_x( 0, std_x );
	normal_distribution< double > dist_y( 0, std_y );
	normal_distribution< double > dist_theta( 0, std_yaw );

	for( int i = 0; i < particles.size(); ++i )
	{
		Particle& p = particles[ i ];

		p.x += dist_x( gen );
		p.y += dist_y( gen );
		p.theta += dist_theta( gen );
	}
#endif
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for( int i = 0; i < observations.size(); ++i )
	{
		LandmarkObs& o = observations[ i ];
		
		double smallestDist = 99999;
		LandmarkObs bestMatch;

		for( int j = 0; j < predicted.size(); ++j )
		{
			const LandmarkObs p = predicted[ j ];

			/* calc Euclidean distance */
			double d = sqrt( pow( o.x - p.x, 2 ) + pow( o.y - p.y, 2 ) );
			if( d < smallestDist )
			{
				bestMatch = p;
				smallestDist = d;
			}
		}
		o = bestMatch;
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


	/* transform map_landmarks (class Map) to landmarks vector to be able to use them in dataAssociation() */
	std::vector< LandmarkObs > vLandmarksMap;
	for( int i = 0; i < map_landmarks.landmark_list.size(); ++i )
	{
		LandmarkObs l = { i, map_landmarks.landmark_list[ i ].x_f, map_landmarks.landmark_list[ i ].y_f };
		vLandmarksMap.push_back( l );
	}

	/* iterate through all particles and update weights */
	for( int i = 0; i < particles.size(); ++i )
	{
		Particle& p = particles[ i ];

		/* transform observations to map coordinates for the current particle */
		std::vector< LandmarkObs > observations_map;
		for( int j = 0; j < observations.size(); ++j )
		{
			const LandmarkObs& o = observations[ j ];

			double x_map = p.x + cos( p.theta ) * o.x - sin( p.theta ) * o.y;
			double y_map = p.y + sin( p.theta ) * o.x - cos( p.theta ) * o.y;

			observations_map.push_back( { j, x_map, y_map } );
		}

		/* associate observations to closest landmarks in the map */
		std::vector< LandmarkObs > observations_map_associations = observations_map;
		dataAssociation( vLandmarksMap, observations_map_associations );

		/* update weight of the particle */
		const double std_x = std_landmark[ 0 ];
		const double std_y = std_landmark[ 1 ];
		const double gauss_norm = ( 1.0 / ( 2 * M_PI * std_x * std_y ) );		// calculate normalization term
		double weight = 1.0f;
		
		for( int j = 0; j < observations_map.size(); ++j )
		{
			const LandmarkObs& o = observations_map[ j ];				// observed landmark in map coordinate system
			const LandmarkObs& u = observations_map_associations[ j ];	// nearest landmark in map coordinate system
			// calculate exponent
			double exponent = pow( ( o.x - u.x ), 2 ) / ( 2 * std_x * std_x ) + pow( ( o.y - u.y ), 2 ) / ( 2 * std_y * std_y );

			//	calculate weight using normalization terms and exponent
			weight = weight * ( gauss_norm * exp( -exponent ) );
		}

		weights[ i ] = weight;
		particles[ i ].weight = weight;
	}

	/* normalizing weights */
	//int max 
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
