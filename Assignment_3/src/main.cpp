////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

// Eigen for matrix operations
#include <Eigen/Dense>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// JSON parser library (https://github.com/nlohmann/json)
#include "json.hpp"
using json = nlohmann::json;

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Define types & classes
////////////////////////////////////////////////////////////////////////////////

struct Ray {
	Vector3d origin;
	Vector3d direction;
	Ray() { }
	Ray(Vector3d o, Vector3d d) : origin(o), direction(d) { }
};

struct Light {
	Vector3d position;
	Vector3d intensity;
};

struct Intersection {
	Vector3d position;
	Vector3d normal;
	double ray_param;
};

struct Camera {
	bool is_perspective;
	Vector3d position;
	double field_of_view; // between 0 and PI
	double focal_length;
	double lens_radius; // for depth of field
};

struct Material {
	Vector3d ambient_color;
	Vector3d diffuse_color;
	Vector3d specular_color;
	double specular_exponent; // Also called "shininess"

	Vector3d reflection_color;
	Vector3d refraction_color;
	double refraction_index;
};

struct Object {
	Material material;
	virtual ~Object() = default; // Classes with virtual methods should have a virtual destructor!
	virtual bool intersect(const Ray &ray, Intersection &hit) = 0;
};

// We use smart pointers to hold objects as this is a virtual class
typedef std::shared_ptr<Object> ObjectPtr;

struct Sphere : public Object {
	Vector3d position;
	double radius;

	virtual ~Sphere() = default;
	virtual bool intersect(const Ray &ray, Intersection &hit) override;
};

struct Parallelogram : public Object {
	Vector3d origin;
	Vector3d u;
	Vector3d v;

	virtual ~Parallelogram() = default;
	virtual bool intersect(const Ray &ray, Intersection &hit) override;
};

struct Scene {
	Vector3d background_color;
	Vector3d ambient_light;

	Camera camera;
	std::vector<Material> materials;
	std::vector<Light> lights;
	std::vector<ObjectPtr> objects;
};

////////////////////////////////////////////////////////////////////////////////

bool Sphere::intersect(const Ray &ray, Intersection &hit) {
	double t = (this->position - ray.origin).dot(ray.direction.normalized());
	Vector3d intersection = ray.origin + t * ray.direction.normalized();
	double r = (intersection - this->position).norm();
	if ( r < this->radius) {
		double d = std::sqrt( this->radius * this->radius
							- (intersection - this->position).squaredNorm());
		hit.position = intersection - d * ray.direction.normalized(); 
		hit.normal = (hit.position - this->position).normalized();
		hit.ray_param = t - d ; //  Not sure what it is, assume it's the scalar.
		return true;
	} else {
		return false;
	}
}

bool Parallelogram::intersect(const Ray &ray, Intersection &hit) {
	// compute intersect coefficient
	Vector3d plain_norm = u.cross(v);
	double t = (this->origin - ray.origin).dot(plain_norm)
		     / (ray.direction(0) + ray.direction(1) + ray.direction(2));
	Vector3d intersection = ray.origin + ray.direction * t;
	MatrixXd p_gram;
	p_gram << u, v;
	Vector2d portion = p_gram.colPivHouseholderQr().solve(intersection - this->origin);
	if (0 < portion(0) && portion(0) < 1 && 0 < portion(1) && portion(1) < 1) {
		hit.normal = u.cross(v).normalized();
		hit.position = intersection;
		hit.ray_param = t;
		return true;
	} else {
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Define ray-tracing functions
////////////////////////////////////////////////////////////////////////////////

// Function declaration here (could be put in a header file)
Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &object, const Intersection &hit, int max_bounce);
Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit);
bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light);
Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce);

// -----------------------------------------------------------------------------

Vector3d ray_color(const Scene &scene, const Ray &ray, const Object &obj, const Intersection &hit, int max_bounce) {
	// Material for hit object
	const Material &mat = obj.material;

	// Ambient light contribution
	Vector3d ambient_color = obj.material.ambient_color.array() * scene.ambient_light.array();

	// Punctual lights contribution (direct lighting)
	Vector3d lights_color(0, 0, 0);
	for (const Light &light : scene.lights) {
		Vector3d Li = (light.position - hit.position).normalized();
		Vector3d N = hit.normal;

		Ray shadow;
		shadow.origin = hit.position;
		shadow.direction = Li;

		if (!is_light_visible(scene, shadow, light)) {
			continue; // skip if in shadow
		}
			// Diffuse contribution
		Vector3d diffuse = mat.diffuse_color * std::max(Li.dot(N), 0.0);

		Vector3d specular(0, 0, 0);
		Vector3d hit_pos = ( Li - ray.direction.normalized()).normalized();
		double phong_exp =  hit_pos.dot(hit.normal); 
		specular = mat.specular_color * std::max(pow(phong_exp, mat.specular_exponent), 0.0);

		// Attenuate lights according to the squared distance to the lights
		Vector3d D = light.position - hit.position;
		lights_color += (diffuse + specular).cwiseProduct(light.intensity) /  D.squaredNorm();
	}
	// TODO: Compute the color of the reflected ray and add its contribution to the current point color.
	Vector3d reflection_color(0, 0, 0);
    if (max_bounce > 0) {
        Ray reflection;
        reflection.origin = hit.position;
        reflection.direction =
            ray.direction.normalized() -
            2 * (hit.normal.dot(ray.direction.normalized())) *
                hit.normal; // compute reflect ray direction
        reflection_color = mat.reflection_color.cwiseProduct(
            shoot_ray(scene, reflection, max_bounce-1));
    }	// TODO: Compute the color of the refracted ray and add its contribution to the current point color.
	//       Make sure to check for total internal reflection before shooting a new ray.
	Vector3d refraction_color(0, 0, 0);
	Ray refraction;
	refraction.origin = hit.position;
	Vector3d ray_direct = ray.direction.normalized();
	double c1 = ray_direct.dot(hit.normal);
	double c2 = sqrt(1-(1/mat.refraction_index)*(1/mat.refraction_index)*(1-(hit.normal.dot(ray_direct)*(hit.normal.dot(ray_direct)))));
	refraction.direction = 1/mat.refraction_index*(ray_direct + c1*hit.normal) - hit.normal*c2;
	Intersection rf_hit;
	Object *r_obj = find_nearest_object(scene, refraction, rf_hit);
	if(r_obj) {
		double c3 = refraction.direction.normalized().dot(rf_hit.normal);
		double c4 = sqrt(1-(mat.refraction_index)*(mat.refraction_index)*(1-(rf_hit.normal.dot(refraction.direction)*(rf_hit.normal.dot(refraction.direction)))));
		Ray ref_ray;
		ref_ray.origin = rf_hit.position;
		ref_ray.direction = mat.refraction_index*(refraction.direction + c3*rf_hit.normal) - rf_hit.normal*c4;
		refraction_color = mat.refraction_color.cwiseProduct(ray_color(scene, refraction, *r_obj, rf_hit, max_bounce - 1));
	}


	// Rendering equation
	Vector3d C = ambient_color + lights_color + reflection_color + refraction_color;
	// Vector3d C = ambient_color + lights_color; 

	return C;
}

// -----------------------------------------------------------------------------

Object * find_nearest_object(const Scene &scene, const Ray &ray, Intersection &closest_hit) {
	int closest_index = -1;
	Intersection hit;
	const double eps = 1e-10; // small epsilon for numerical precision
	double closest_param = 100; // a large number that is greater than t
	long int count = 0;
	for (const ObjectPtr &obj: scene.objects) {
		if (obj->intersect(ray, hit)) {
			if (hit.ray_param < closest_param && hit.ray_param > eps) {
				closest_index = count;
				closest_hit = hit;
				closest_param = hit.ray_param; // the distance to ray origin
			}
		}
		count++;
	}
	if (closest_index < 0) {
		// Return a NULL pointer
		return nullptr;
	} else {
		// Return a pointer to the hit object. Don't forget to set 'closest_hit' accordingly!
		return scene.objects[closest_index].get();
	}
}

bool is_light_visible(const Scene &scene, const Ray &ray, const Light &light) {
	// TODO: Determine if the light is visible here
	const double eps = 1e-10; // small epsilon for numerical precision
	Intersection hit_obj;
	double light_param = (light.position(0) - ray.origin(0)) / ray.direction(0);
	for (const ObjectPtr &h_obj : scene.objects) {
		if (h_obj->intersect(ray, hit_obj) && hit_obj.ray_param > eps && hit_obj.ray_param < light_param) {
			return false;
		}
	}
	return true;
}

Vector3d shoot_ray(const Scene &scene, const Ray &ray, int max_bounce) {
	Intersection hit;
	if (Object *obj = find_nearest_object(scene, ray, hit)) {
		// 'obj' is not null and points to the object of the scene hit by the ray
		return ray_color(scene, ray, *obj, hit, max_bounce);
	} else {
		// 'obj' is null, we must return the background color
		return scene.background_color;
	}
}

////////////////////////////////////////////////////////////////////////////////

void render_scene(const Scene &scene) {
	std::cout << "Simple ray tracer." << std::endl;

	int w = 640;
	int h = 480;
	MatrixXd R = MatrixXd::Zero(w, h);
	MatrixXd G = MatrixXd::Zero(w, h);
	MatrixXd B = MatrixXd::Zero(w, h);
	MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

	// The camera always points in the direction -z
	// The sensor grid is at a distance 'focal_length' from the camera center,
	// and covers an viewing angle given by 'field_of_view'.
	double aspect_ratio = double(w) / double(h);
	double scale_y = scene.camera.focal_length * tan(scene.camera.field_of_view/2) * 2;
	double scale_x = scale_y * aspect_ratio;//

	// The pixel grid through which we shoot rays is at a distance 'focal_length'
	// from the sensor, and is scaled from the canonical [-1,1] in order
	// to produce the target field of view.
	int iter = 150;
	Vector3d grid_origin(-scale_x, scale_y, -scene.camera.focal_length);
	Vector3d x_displacement(2.0/w*scale_x, 0, 0);
	Vector3d y_displacement(0, -2.0/h*scale_y, 0);

	for (unsigned i = 0; i < w; ++i) {
		for (unsigned j = 0; j < h; ++j) {
			for (unsigned k = 0; k < iter; ++k) {
				Vector3d shift = grid_origin + (i+0.5)*x_displacement + (j+0.5)*y_displacement;

				// Prepare the ray
				Ray ray;
				double origin_x_offset = (2 *(float)rand() / RAND_MAX - 1);
				double origin_y_offset = (2 *(float)rand() / RAND_MAX - 1) * sqrt( 1 - origin_x_offset * origin_x_offset);
				Vector3d offset(origin_x_offset, origin_y_offset,0);

				if (scene.camera.is_perspective) {
					// Perspective camera
					ray.origin = scene.camera.position + offset * scene.camera.lens_radius;
					ray.direction = shift - ray.origin;
				} else {
					// Orthographic camera
					ray.origin = scene.camera.position + Vector3d(shift[0], shift[1], 0);
					ray.direction = Vector3d(0, 0, -1);
				}

				int max_bounce = 5;
				Vector3d C = shoot_ray(scene, ray, max_bounce);
				R(i, j) += C(0);
				G(i, j) += C(1);
				B(i, j) += C(2);
			}
				A(i, j) = 1;
		}
	}

	// Save to png
	const std::string filename("raytrace.png");
	write_matrix_to_png(R/iter, G/iter, B/iter, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

Scene load_scene(const std::string &filename) {
	Scene scene;

	// Load json data from scene file
	json data;
	std::ifstream in(filename);
	in >> data;

	// Helper function to read a Vector3d from a json array
	auto read_vec3 = [] (const json &x) {
		return Vector3d(x[0], x[1], x[2]);
	};

	// Read scene info
	scene.background_color = read_vec3(data["Scene"]["Background"]);
	scene.ambient_light = read_vec3(data["Scene"]["Ambient"]);

	// Read camera info
	scene.camera.is_perspective = data["Camera"]["IsPerspective"];
	scene.camera.position = read_vec3(data["Camera"]["Position"]);
	scene.camera.field_of_view = data["Camera"]["FieldOfView"];
	scene.camera.focal_length = data["Camera"]["FocalLength"];
	scene.camera.lens_radius = data["Camera"]["LensRadius"];

	// Read materials
	for (const auto &entry : data["Materials"]) {
		Material mat;
		mat.ambient_color = read_vec3(entry["Ambient"]);
		mat.diffuse_color = read_vec3(entry["Diffuse"]);
		mat.specular_color = read_vec3(entry["Specular"]);
		mat.reflection_color = read_vec3(entry["Mirror"]);
		mat.refraction_color = read_vec3(entry["Refraction"]);
		mat.refraction_index = entry["RefractionIndex"];
		mat.specular_exponent = entry["Shininess"];
		scene.materials.push_back(mat);
	}

	// Read lights
	for (const auto &entry : data["Lights"]) {
		Light light;
		light.position = read_vec3(entry["Position"]);
		light.intensity = read_vec3(entry["Color"]);
		scene.lights.push_back(light);
	}

	// Read objects
	for (const auto &entry : data["Objects"]) {
		ObjectPtr object;
		if (entry["Type"] == "Sphere") {
			auto sphere = std::make_shared<Sphere>();
			sphere->position = read_vec3(entry["Position"]);
			sphere->radius = entry["Radius"];
			object = sphere;
		} else if (entry["Type"] == "Parallelogram") {
			// TODO
		}
		object->material = scene.materials[entry["Material"]];
		scene.objects.push_back(object);
	}

	return scene;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " scene.json" << std::endl;
		return 1;
	}
	Scene scene = load_scene(argv[1]);
	render_scene(scene);
	return 0;
}
