////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <queue>

// Eigen for matrix operations
#include <Eigen/Dense>
#include <Eigen/Geometry>

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

struct AABBTree {
	struct Node {
		AlignedBox3d bbox;
		int parent; // Index of the parent node (-1 for root)
		int left; // Index of the left child (-1 for a leaf)
		int right; // Index of the right child (-1 for a leaf)
		int triangle; // Index of the node triangle (-1 for internal nodes)
	};

	std::vector<Node> nodes;
	int root;

	AABBTree() = default; // Default empty constructor
	AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh
};

struct Mesh : public Object {
	MatrixXd vertices; // n x 3 matrix (n points)
	MatrixXi facets; // m x 3 matrix (m triangles)

	AABBTree bvh;

	Mesh() = default; // Default empty constructor
	Mesh(const std::string &filename);
	virtual ~Mesh() = default;
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

// Read a triangle mesh from an off file
void load_off(const std::string &filename, MatrixXd &V, MatrixXi &F) {
	std::ifstream in(filename);
	std::string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	V.resize(nv, 3);
	F.resize(nf, 3);
	for (int i = 0; i < nv; ++i) {
		in >> V(i, 0) >> V(i, 1) >> V(i, 2);
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> F(i, 0) >> F(i, 1) >> F(i, 2);
		assert(s == 3);
	}
}

Mesh::Mesh(const std::string &filename) {
	// Load a mesh from a file (assuming this is a .off file), and create a bvh
	load_off(filename, vertices, facets);
	bvh = AABBTree(vertices, facets);
}

////////////////////////////////////////////////////////////////////////////////
// BVH Implementation
////////////////////////////////////////////////////////////////////////////////

// Bounding box of a triangle
AlignedBox3d bbox_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c) {
	AlignedBox3d box;
	box.extend(a);
	box.extend(b);
	box.extend(c);
	return box;
}

// centroid struct for sorting
struct centroid {
	double x; 
	double y; 
	double z; 
	long idx; // index
} centroid;

// sorting criteria.
bool sort_x(struct centroid* left, struct centroid* right) {
			return left->x < right->x;
}

bool sort_y(struct centroid* left, struct centroid* right) {
			return left->y < right->y;
}

bool sort_z(struct centroid* left, struct centroid* right) {
			return left->z < right->z;
}

AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F) {
	// Compute the centroids of all the triangles in the input mesh
	MatrixXd centroids(F.rows(), V.cols());
	centroids.setZero();
	for (int i = 0; i < F.rows(); ++i) {
		for (int k = 0; k < F.cols(); ++k) {
			centroids.row(i) += V.row(F(i, k));
		}
		centroids.row(i) /= F.cols();
	}
	// create vector, associate centroid with index
	long curr = 0; // for parent index
	std::queue<std::vector<struct centroid*>*> q_centroid;
	std::queue<int> q_parent;
	q_parent.push(-1);
	// create vector, associate centroid with index
	std::vector<struct centroid*>* v_centroids = new std::vector<struct centroid*>;
	q_centroid.push(v_centroids);
	for (int i = 0; i < centroids.rows(); ++i) {
		struct centroid* s_centroid = new struct centroid;
		s_centroid->x = centroids(i, 0);
		s_centroid->y = centroids(i, 1);
		s_centroid->z = centroids(i, 2);
		s_centroid->idx = i;
		v_centroids->push_back(s_centroid);
	}
	while (q_centroid.size() > 0) {
		std::vector<struct centroid*>* list = q_centroid.front();
		int parent = q_parent.front();
		q_centroid.pop();
		q_parent.pop();
		if (list->size() == 1) {
			struct centroid* triangle = list->front();
			list->pop_back();
			Node *node = new Node;
			node->left = -1;
			node->right = -1;
			node->triangle = triangle->idx;
			node->parent = parent;
			nodes.push_back(*node);
			++curr; // move index forward
		} else {
			switch (curr % 3) { // change sorting direction
				case 0:
					std::sort(list->begin(), list->end(), sort_x);
					break;
				case 1:
					std::sort(list->begin(), list->end(), sort_y);
					break;
				case 2:
					std::sort(list->begin(), list->end(), sort_z);
					break;
			}
			int mid = list->size() / 2;
			std::vector<struct centroid*>* left = new std::vector<struct centroid*>;
			std::vector<struct centroid*>* right = new std::vector<struct centroid*>;
			for (int i = 0; i< mid; ++i) {
				left->push_back((*list)[i]);
			}
			for (int i = mid; i < list->size(); ++i) {
				right->push_back((*list)[i]);
			}
			Node* node = new Node;
			node->parent = parent;
			node->triangle = -1;
			node->left = -1;
			node->right = -1;
			q_centroid.push(left);
			q_centroid.push(right);
			q_parent.push(curr);
			q_parent.push(curr);
			nodes.push_back(*node);
			++curr; // move index forward
		}
	}
	// build left & right
	for ( int i = nodes.size() - 1; i > 0; --i) {
		if (nodes.at(i).triangle != -1) { // Build bbox
			nodes.at(i).bbox = bbox_triangle(V.row(F(nodes.at(i).triangle, 0)),
				V.row(F(nodes.at(i).triangle, 1)), V.row(F(nodes.at(i).triangle, 2)));
		} else {
			nodes.at(i).bbox = nodes.at(nodes.at(i).left).bbox.merged(nodes.at(nodes.at(i).left).bbox);
		}
		int parent = nodes.at(i).parent;
		if (nodes.at(parent).left == -1) {
			nodes.at(parent).left = i;
		} else {
			nodes.at(parent).right = i;
		}
	}
	// for test
	// for (std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
		// printf("%d, %d, %d\n ", it->parent, it->left, it->right);
	// }
	// delete centroids
}

////////////////////////////////////////////////////////////////////////////////

bool Sphere::intersect(const Ray &ray, Intersection &hit) {
	// TODO (Assignment 2)
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
	return false;
}

bool Parallelogram::intersect(const Ray &ray, Intersection &hit) {
	// TODO (Assignment 2)
	Vector3d plain_norm = u.cross(v);
	double t = (this->origin - ray.origin).dot(plain_norm)
		     / (ray.direction.dot(plain_norm));
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

// -----------------------------------------------------------------------------

bool intersect_triangle(const Ray &ray, const Vector3d &a, const Vector3d &b, const Vector3d &c, Intersection &hit) {
	// TODO (Assignment 3)
	//
	// Compute whether the ray intersects the given triangle.
	// If you have done the parallelogram case, this should be very similar to it.
	Vector3d u = b - a;
	Vector3d v = c - a;
	Vector3d plain_norm = u.cross(v).normalized();
	double t = (a - ray.origin).dot(plain_norm) / ray.direction.dot(plain_norm);
	Vector3d intersection = ray.origin + ray.direction * t;
	MatrixXd triangle(3, 2);
	triangle.col(0) = u; 
	triangle.col(1) = v; 
	Vector2d portion = triangle.colPivHouseholderQr().solve(intersection - a);
	if (0 < portion(0) &&  0 < portion(1) && portion(0) + portion(1) < 1) {
		hit.normal = u.cross(v).normalized();
		hit.position = intersection;
		hit.ray_param = t;
		return true;
	} else {
		return false;
	}
}

bool intersect_box(const Ray &ray, const AlignedBox3d &box) {
	// TODO (Assignment 3)
	//
	// Compute whether the ray intersects the given box.
	// There is no need to set the resulting normal and ray parameter, since
	// we are not testing with the real surface here anyway.
	double scalar_front = (box.corner(box.BottomLeft)(1) - ray.origin(1)) / ray.direction(1); 
	double scalar_back = (box.corner(box.TopLeft)(1) - ray.origin(1)) / ray.direction(1); 
	Vector3d intersect_1 = ray.origin + scalar_front * ray.direction;
	Vector3d intersect_2 = ray.origin + scalar_back * ray.direction;
	if ((box.corner(AlignedBox3d::TopLeft)(0) <= intersect_1(0) && box.corner(box.BottomRightCeil)(0) >= intersect_1(0)
	&& box.corner(box.TopLeft)(1) <= intersect_1(1) && box.corner(box.BottomRightCeil)(1) >= intersect_1(1)
	&& box.corner(box.TopLeft)(2) <= intersect_1(2) && box.corner(box.BottomRightCeil)(2) >= intersect_1(2)) 
	|| (box.corner(box.TopLeft)(0) <= intersect_2(0) && box.corner(box.BottomRightCeil)(0) >= intersect_2(0)
	&& box.corner(box.TopLeft)(1) <= intersect_2(1) && box.corner(box.BottomRightCeil)(1) >= intersect_2(1)
	&& box.corner(box.TopLeft)(2) <= intersect_2(2) && box.corner(box.BottomRightCeil)(2) >= intersect_2(2))) {
		return true;
	}
	return false;
}

bool Mesh::intersect(const Ray &ray, Intersection &closest_hit) {
	// TODO (Assignment 3)

	// Method (1): Traverse every triangle and return the closest hit.
	// Method (1):
	Intersection hit;
	double min_param = 1000;
	bool intersect = false;
	for (int i = 0; i < facets.rows(); ++i){
		int a = facets(i, 0);
		int b = facets(i, 1);
		int c = facets(i, 2);
		Vector3d v_a(vertices(a, 0), vertices(a, 1), vertices(a, 2));
		Vector3d v_b(vertices(b, 0), vertices(b, 1), vertices(b, 2));
		Vector3d v_c(vertices(c, 0), vertices(c, 1), vertices(c, 2));
		if (intersect_triangle(ray, v_a, v_b, v_c, hit)) {
			if (hit.ray_param < min_param) {
				closest_hit = hit;
				min_param = hit.ray_param;
				intersect = true;
			}
		}
	}
	// Method (2): Traverse the BVH tree and test the intersection with a
	// triangles at the leaf nodes that intersects the input ray.
	return intersect;
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

		// TODO (Assignment 2, specular contribution)
		Vector3d specular(0, 0, 0);
		Vector3d hit_pos = ( Li - ray.direction.normalized()).normalized();
		double phong_exp = hit_pos.dot(hit.normal); 
		specular = mat.specular_color * std::max(pow(phong_exp, mat.specular_exponent), 0.0);

		// Attenuate lights according to the squared distance to the lights
		Vector3d D = light.position - hit.position;
		lights_color += (diffuse + specular).cwiseProduct(light.intensity) /  D.squaredNorm();
	}
	// TODO (Assignment 2, reflected ray)
	Vector3d reflection_color(0, 0, 0);
	if (max_bounce > 0) {
		Ray reflection;
		reflection.origin = hit.position;
		reflection.direction = ray.direction.normalized() - 2 * (hit.normal.dot(ray.direction.normalized())) * hit.normal; // compute reflect ray direction
		Intersection r_hit;
		Object *n_obj = find_nearest_object(scene, reflection, r_hit);
		if (n_obj) { // make sure the intersection point is in the positive direction of reflection
			reflection_color = mat.reflection_color.cwiseProduct(ray_color(scene, reflection, *n_obj, r_hit, max_bounce - 1));
		} 
	}

	// TODO (Assignment 2, refracted ray)
	Vector3d refraction_color(0, 0, 0);

	// Rendering equation
	Vector3d C = ambient_color + lights_color + reflection_color + refraction_color;

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
	// TODO (Assignment 2, shadow ray)
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
	if (Object * obj = find_nearest_object(scene, ray, hit)) {
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
	int iter = 1;
	Vector3d grid_origin(-scale_x, scale_y, -scene.camera.focal_length);
	Vector3d x_displacement(2.0/w*scale_x, 0, 0);
	Vector3d y_displacement(0, -2.0/h*scale_y, 0);

	for (unsigned i = 0; i < w; ++i) {
		std::cout << std::fixed << std::setprecision(2);
		std::cout << "Ray tracing: " << (100.0 * i) / w << "%\r" << std::flush;
		for (unsigned j = 0; j < h; ++j) {
			for (unsigned k = 0; k < iter; ++k) {
				Vector3d shift = grid_origin + (i+0.5)*x_displacement + (j+0.5)*y_displacement;

				// Prepare the ray
				Ray ray;
				double origin_x_offset = (2 *(float)rand() / RAND_MAX - 1);
				double origin_y_offset = (2 *(float)rand() / RAND_MAX - 1) * sqrt( 1 - origin_x_offset * origin_x_offset);
				Vector3d offset(origin_x_offset, origin_y_offset,0);
			// TODO (Assignment 2, depth of field)
				// Prepare the ray

				if (scene.camera.is_perspective) {
					// Perspective camera
					// TODO (Assignment 2, perspective camera)
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

	std::cout << "Ray tracing: 100%  " << std::endl;

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
		} else if (entry["Type"] == "Mesh") {
			// Load mesh from a file
			std::string filename = std::string(DATA_DIR) + entry["Path"].get<std::string>();
			object = std::make_shared<Mesh>(filename);
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
