// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"


// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere() {
	std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

	const std::string filename("sphere_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			// NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const double sphere_radius = 0.9;

			if (ray_on_xy.norm() < sphere_radius) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

void raytrace_parallelogram() {
	std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

	const std::string filename("plane_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.6, -0.2, 0);
	Vector3d pgram_u(0.3, -0.2, 0);
	Vector3d pgram_v(-0.3, 0.4, -1);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// TODO: Check if the ray intersects with the parallelogram
			Matrix2d pgram_point;
			pgram_point <<  pgram_u(0) - pgram_origin(0), pgram_v(0) - pgram_origin(0),
							pgram_u(1) - pgram_origin(1), pgram_v(1) - pgram_origin(1);
			Vector2d ray_origin_2d(ray_origin(0) - pgram_origin(0), ray_origin(1) - pgram_origin(1));
			Vector2d intersect = pgram_point.colPivHouseholderQr().solve(ray_origin_2d);
			if (0 <= intersect(0) 
				&& intersect(0) <= 1
				&& 0 <= intersect(1)
				&& intersect(1) <= 1) {
				// TODO: The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection(0,0,0);
				ray_intersection = pgram_origin 
								+ intersect(0) * pgram_u 
								+ intersect(1) * pgram_v;
				// TODO: Compute normal at the intersection point
				
				Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_perspective() {
	std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

	const std::string filename("plane_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d cam_pos(0, 0, 2); // assume camera position.
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-1, -1, -1);
	Vector3d pgram_u(1, -0.4, -1);
	Vector3d pgram_v(-1, 1, -1.2);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			Vector3d pixel = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = pixel - cam_pos;
			Vector3d pgram_v_u = pgram_u - pgram_origin;
			Vector3d pgram_v_v = pgram_v - pgram_origin;
			Vector3d norm = pgram_v_u.cross(pgram_v_v).normalized();
			double t = ( norm(0) * (pgram_origin(0) - cam_pos(0))
					   + norm(1) * (pgram_origin(1) - cam_pos(1))
					   + norm(2) * (pgram_origin(2) - cam_pos(2)) )
					   / (norm(0) * ray_direction(0) 
					   +  norm(1) * ray_direction(1)
					   +  norm(2) * ray_direction(2));  // scalar of the parametrized line
			Vector3d intersect_point = cam_pos + t * ray_direction;
			Vector3d intersect_vector = intersect_point - pgram_origin;
			MatrixXd pgram_vectors(3, 2);
			pgram_vectors << pgram_v_u, pgram_v_v;
			Vector2d intersect = pgram_vectors.colPivHouseholderQr().solve(intersect_vector);
			if (0 <= intersect(0) 
				&& intersect(0) <= 1
				&& 0 <= intersect(1)
				&& intersect(1) <= 1) {
				Vector3d ray_intersection(0,0,0);
				ray_intersection = intersect_point;
				Vector3d ray_normal = pgram_v_u.cross(pgram_v_v).normalized();
				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);
	double ambient = 0.1;
	MatrixXd diffuse = MatrixXd::Zero(800, 800);
	MatrixXd specular = MatrixXd::Zero(800, 800);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			// NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const double sphere_radius = 0.9;

			if (ray_on_xy.norm() < sphere_radius) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// TODO: Add shading parameter here
				diffuse(i,j) = std::max((light_position-ray_intersection).normalized().dot(ray_normal), (double)0);
				// Compute Specular shading
				Vector3d h = ((light_position - ray_intersection).normalized() - ray_direction).normalized();
				double phong_exp = h.transpose() * ray_normal;
				phong_exp = std::max(phong_exp, (double)0);
				// Simple diffuse model
				C(i,j) = ambient + diffuse(i,j) + pow(phong_exp, 1000);

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C/2.5,C/1.25,C, A, filename);
}

int main() {
	raytrace_sphere();
	raytrace_parallelogram();
	raytrace_perspective();
	raytrace_shading();

	return 0;
}
