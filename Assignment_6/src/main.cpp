// This example is heavily based on the tutorial at https://open.gl

////////////////////////////////////////////////////////////////////////////////
// OpenGL Helpers to reduce the clutter
#include "helpers.h"
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
// Linear Algebra Library
#include <Eigen/Dense>
#include <Eigen/Geometry>
// STL headers
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

int obj_num = 0; // track obj number


enum MODE {
	WIREFRAME = 0,
	FLAT,
	PHONG,
};

void cube_init(Eigen::MatrixXf& cube_V) {
	// init Vertex
	cube_V.resize(36, 3);
	cube_V <<
		-0.5, -0.5, -0.5,
		-0.5, -0.5, 0.5,
		0.5, -0.5, 0.5,

		-0.5, -0.5, -0.5,
		0.5, -0.5, 0.5,
		0.5, -0.5, -0.5,

		-0.5, -0.5, -0.5,
		0.5, 0.5, -0.5,
		-0.5, 0.5, -0.5,

		-0.5, -0.5, -0.5,
		0.5, -0.5, -0.5,
		0.5, 0.5, -0.5,

		-0.5, -0.5, -0.5,
		-0.5, 0.5, -0.5,
		-0.5, 0.5, 0.5,

		-0.5, -0.5, -0.5,
		-0.5, 0.5, 0.5,
		-0.5, -0.5, 0.5,

		0.5, 0.5, 0.5,
		0.5, -0.5, -0.5,
		0.5, -0.5, 0.5,

		0.5, 0.5, 0.5,
		0.5, 0.5, -0.5,
		0.5, -0.5, -0.5,

		0.5, 0.5, 0.5,
		0.5, -0.5, 0.5,
		-0.5, -0.5, 0.5,

		0.5, 0.5, 0.5,
		-0.5, -0.5, 0.5,
		-0.5, 0.5, 0.5,

		0.5, 0.5, 0.5,
		-0.5, 0.5, 0.5,
		-0.5, 0.5, -0.5,

		0.5, 0.5, 0.5,
		-0.5, 0.5, -0.5,
		0.5, 0.5, -0.5;

	cube_V.transposeInPlace();
}


////////////////////////////////////////////////////////////////////////////////


// Mesh object, with both CPU data (Eigen::Matrix) and GPU data (the VBOs)
struct Mesh {
	Eigen::MatrixXf V; // mesh vertices [3 x n]

	Eigen::MatrixXf N; // mesh trasformation matrix [4 x j]
	Eigen::MatrixXf NV; // mesh transformation index
	std::vector<Eigen::Matrix4f*> mtx_list; // track transform mtx
	std::vector<int> obj_v_num;  // track how many triangle per obj 
	std::vector<int> off_num = { 0 };  // track offset
	std::vector<int> obj_select; // track whether selected
	MODE mode = WIREFRAME; // display mode

	// VBO storing vertex position attributes
	VertexBufferObject V_vbo;

	// VBO storing normals
	VertexBufferObject N_vbo;

	// VBO storing normals of vertices
	VertexBufferObject NV_vbo;

	// VAO storing the layout of the shader program for the object 'bunny'
	VertexArrayObject vao;
};

Mesh bunny;

////////////////////////////////////////////////////////////////////////////////

// Read a triangle mesh from an off file
void load_off(const std::string& filename, Eigen::MatrixXf& V) {
	std::ifstream in(filename);
	std::string token;
	Eigen::MatrixXf tmpV;
	Eigen::MatrixXi tmpF;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	tmpV.resize(3, nv);
	tmpF.resize(3, nf);
	V.resize(3, nf * 3);

	for (int i = 0; i < nv; ++i) {
		in >> tmpV(0, i) >> tmpV(1, i) >> tmpV(2, i);
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> tmpF(0, i) >> tmpF(1, i) >> tmpF(2, i);
		assert(s == 3);
	}

	// add normal


	for (int i = 0; i < nf; ++i) {
		V.col(3 * i) = tmpV.col(tmpF(0, i));
		V.col(3 * i + 1) = tmpV.col(tmpF(1, i));
		V.col(3 * i + 2) = tmpV.col(tmpF(2, i));
	}

	long orig_N_size = bunny.N.cols();
	bunny.N.conservativeResize(bunny.N.rows(), bunny.N.cols() + 3 * nf);
	for (int i = 0; i < 3 * nf; i += 3) {
		Eigen::Vector3f a = V.col(i + 1) - V.col(i);
		Eigen::Vector3f b = V.col(i + 2) - V.col(i);
		bunny.N.col(orig_N_size + i) = bunny.N.col(orig_N_size + i + 1) = bunny.N.col(orig_N_size + i + 2)
			= a.cross(b).normalized();
	}
	bunny.N_vbo.update(bunny.N);

	// compute vector normal
	std::vector<Eigen::Vector3f> tmpNV;
	for (int i = 0; i < nv; ++i) {
		tmpNV.push_back(Eigen::Vector3f(0, 0, 0));
		for (int j = 0; j < nf; ++j) {
			if (tmpF(0, j) == i) {
				tmpNV[i] += bunny.N.col(orig_N_size + 3 * j);
			}
			if (tmpF(1, j) == i) {
				tmpNV[i] += bunny.N.col(orig_N_size + 3 * j + 1);
			}
			if (tmpF(2, j) == i) {
				tmpNV[i] += bunny.N.col(orig_N_size + 3 * j + 2);
			}
		}
		tmpNV[i].normalize();
	}

	// add to NV list
	long orig_NV_size = bunny.NV.cols();
	bunny.NV.conservativeResize(bunny.NV.rows(), bunny.NV.cols() + 3 * nf);
	for (int i = 0; i < nf; ++i) {
		bunny.NV.col(orig_NV_size + 3 * i) = tmpNV[tmpF(0, i)];
		bunny.NV.col(orig_NV_size + 3 * i + 1) = tmpNV[tmpF(1, i)];
		bunny.NV.col(orig_NV_size + 3 * i + 2) = tmpNV[tmpF(2, i)];
	}
	bunny.NV_vbo.update(bunny.NV);
}

////////////////////////////////////////////////////////////////////////////////

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	// Get viewport size (canvas in number of pixels)
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);

	// Get the size of the window (may be different than the canvas size on retina displays)
	int width_window, height_window;
	glfwGetWindowSize(window, &width_window, &height_window);

	// Get the position of the mouse in the window
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	// Deduce position of the mouse in the viewport
	double highdpi = (double)width / (double)width_window;
	xpos *= highdpi;
	ypos *= highdpi;

	// Convert screen position to the canonical viewing volume
	double xcanonical = ((xpos / double(width)) * 2) - 1;
	double ycanonical = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	// TODO: Ray-casting for object selection (Ex.3)
}

void add_off(std::string dir) {
	Eigen::MatrixXf Tmp_V;
	std::string full_dir = DATA_DIR + dir;
	load_off(full_dir, Tmp_V);

	// add offset
	long orig_V_cols = bunny.V.cols();
	bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + Tmp_V.cols()); // resize for cube
	for (int j = 0; j < Tmp_V.cols(); ++j) {
		bunny.V.col(orig_V_cols + j) = Tmp_V.col(j);
	}
	bunny.V_vbo.update(bunny.V);

	// add model
	float scale = 0;
	float x_off = -(Tmp_V.row(0).minCoeff() + Tmp_V.row(0).maxCoeff()) / 2;
	float y_off = -(Tmp_V.row(1).minCoeff() + Tmp_V.row(1).maxCoeff()) / 2;
	float z_off = -(Tmp_V.row(2).minCoeff() + Tmp_V.row(2).maxCoeff()) / 2;
	float x_len = Tmp_V.row(0).maxCoeff() - Tmp_V.row(0).minCoeff();
	float y_len = Tmp_V.row(1).maxCoeff() - Tmp_V.row(1).minCoeff();
	float z_len = Tmp_V.row(2).maxCoeff() - Tmp_V.row(2).minCoeff();
	if (scale < x_len) {
		scale = x_len;
	}
	if (scale < y_len) {
		scale = y_len;
	}
	if (scale < z_len) {
		scale = z_len;
	}
	Eigen::Matrix4f* trans = new Eigen::Matrix4f; // transformation matrix
	*trans = Eigen::Matrix4f::Identity(); // add identity
	trans->col(3) << x_off, y_off, z_off, 1;
	Eigen::Matrix4f scale_mtx = Eigen::Matrix4f::Identity() / scale; // add identity
	scale_mtx(3, 3) = 1;
	*trans = scale_mtx * (*trans);
	bunny.mtx_list.push_back(trans);
	bunny.obj_v_num.push_back(Tmp_V.cols());
	bunny.off_num.push_back(bunny.off_num.back() + Tmp_V.cols());
	bunny.obj_select.push_back(0); // not selected
	obj_num++; // increment object count
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS) {
		switch (key) {

		case GLFW_KEY_1: {// add cube
			Eigen::MatrixXf cube_V;
			cube_init(cube_V);

			// add normal
			for (int i = 0; i < cube_V.cols(); i += 3) {
				Eigen::Vector3f a = cube_V.col(i + 1) - cube_V.col(i);
				Eigen::Vector3f b = cube_V.col(i + 2) - cube_V.col(i + 1);
				long orig_col = bunny.N.cols();
				bunny.N.conservativeResize(bunny.N.rows(), bunny.N.cols() + 3);
				bunny.N.col(orig_col) = bunny.N.col(orig_col + 1) = bunny.N.col(orig_col + 2) = -1 * a.cross(b).normalized();
			}
			bunny.N_vbo.update(bunny.N);

			// add vector normal

			long orig_NV_cols = bunny.NV.cols();
			long orig_V_cols = bunny.V.cols();
			bunny.NV.conservativeResize(bunny.NV.rows(), bunny.NV.cols() + cube_V.cols()); // resize for cube
			bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + cube_V.cols()); // resize for cube
			for (int j = 0; j < cube_V.cols(); ++j) {
				bunny.V.col(orig_V_cols + j) = cube_V.col(j);
				bunny.NV.col(orig_NV_cols + j) = cube_V.col(j).normalized();
			}
			bunny.NV_vbo.update(bunny.NV);
			bunny.V_vbo.update(bunny.V);

			Eigen::Matrix4f* trans = new Eigen::Matrix4f; // transformation matrix
			*trans = Eigen::Matrix4f::Identity(); // add model 
			bunny.mtx_list.push_back(trans);
			bunny.obj_v_num.push_back(cube_V.cols());
			bunny.off_num.push_back(bunny.off_num.back() + cube_V.cols());
			bunny.obj_select.push_back(0); // not selected
			obj_num++; // increment object count
			break;
		}

		case GLFW_KEY_2: {
			std::string dir = "bumpy_cube.off";
			add_off(dir);
			break;
		}


		case GLFW_KEY_3: {
			std::string dir = "bunny.off";
			add_off(dir);
			break;
		}

		case GLFW_KEY_4: {
			bunny.V.resize(3, 0);
			bunny.V_vbo.update(bunny.V);
			for (int i = 0; i < obj_num; ++i) {
				delete bunny.mtx_list[i];
			}
			bunny.mtx_list.clear();
			bunny.obj_v_num.clear();
			bunny.off_num = { 0 };
			obj_num = 0; // reset 
			break;
		}

		case GLFW_KEY_Q: {

			bunny.mode = WIREFRAME;
			break;
		}

		case GLFW_KEY_W: {

			bunny.mode = FLAT;
			break;
		}

		case GLFW_KEY_E: {

			bunny.mode = PHONG;
			break;
		}

		case GLFW_KEY_P: {
			printf("rows: %ld\n", bunny.V.cols());
			printf("obj_count: %d\n", obj_num);
			std::cout << "Offset:" << std::endl;
			for (int i = 0; i < bunny.off_num.size(); ++i)
				std::cout << bunny.off_num[i] << " ";
			std::cout << std::endl;
			std::cout << "Vertex num:" << std::endl;
			for (int i = 0; i < bunny.obj_v_num.size(); ++i)
				std::cout << bunny.obj_v_num[i] << " ";
			std::cout << std::endl;
			for (int i = 0; i < bunny.mtx_list.size(); ++i)
				std::cout << *(bunny.mtx_list[i]) << "\n" << std::endl;
			printf("mode %d\n", bunny.mode);
			std::cout << bunny.NV << std::endl;
			break;
		}
		default:
			break;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

int main(void) {
	// Initialize the GLFW library
	if (!glfwInit()) {
		return -1;
	}

	// Activate supersampling
	glfwWindowHint(GLFW_SAMPLES, 8);

	// Ensure that we get at least a 3.2 context
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

	// On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

	// Create a windowed mode window and its OpenGL context
	GLFWwindow* window = glfwCreateWindow(640, 640, "[Float] Hello World", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	// Make the window's context current
	glfwMakeContextCurrent(window);

	// Load OpenGL and its extensions
	if (!gladLoadGL()) {
		printf("Failed to load OpenGL and its extensions");
		return(-1);
	}
	printf("OpenGL Version %d.%d loaded", GLVersion.major, GLVersion.minor);

	int major, minor, rev;
	major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
	minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
	rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
	printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
	printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
	printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader = R"(
		#version 150 core

		uniform mat4 model;
		uniform mat4 view;
		uniform mat4 proj;

		in vec3 position;
		in vec3 N;
		in vec3 NV;

		out vec3 N_v;
		out vec3 NV_v;
		out vec3 frag_pos;

		void main() {
			N_v  = normalize(vec3(transpose(inverse(model)) * vec4(N, 1))); 
			NV_v  = normalize(vec3(transpose(inverse(model)) * vec4(NV, 1))); 
			gl_Position = proj * view * model * vec4(position, 1.0);
			frag_pos = vec3(model * vec4(position, 1.0));
		}
	)";

	const GLchar* fragment_shader = R"(
		#version 150 core

		uniform vec3 triangleColor;
		uniform vec3 light_pos;
		uniform int mode;

		in vec3 N_v;
		in vec3 NV_v;
		in vec3 frag_pos;
		out vec4 outColor;

		void main() {
			vec3 light_dir = normalize(frag_pos - light_pos);
			if (mode == 0) {
				outColor = vec4(triangleColor, 1.0);
			} else if (mode == 1) {
				// use N_v for FLAT SHADING
				float diff_mag = max(dot(N_v, light_dir), 0.0);
				vec3 ambient = 0.1 * triangleColor;
				vec3 diffuse = diff_mag * triangleColor;	
				outColor = vec4(ambient + diffuse, 1.0);
			} else if (mode == 2) {
				// use NV_v for PHONG SHADING
				float diff_mag = max(dot(NV_v, light_dir), 0.0);
				// vec3 ambient = 0.1 * triangleColor;
				vec3 diffuse = diff_mag * triangleColor;	
				outColor = vec4(diffuse, 1.0);
			} else {
				// should not reach here
				outColor = vec4(1.0, 1.0, 1.0, 1.0);
			}
		}
	)";

	// Compile the two shaders and upload the binary to the GPU
	// Note that we have to explicitly specify that the output "slot" called outColor
	// is the one that we want in the fragment buffer (and thus on screen)
	program.init(vertex_shader, fragment_shader, "outColor");

	// Prepare a dummy bunny object
	// We need to initialize and fill the two VBO (vertex positions + indices),
	// and use a VAO to store their layout when we use our shader program later.
	{
		// Initialize the VBOs
		bunny.V_vbo.init(GL_FLOAT, GL_ARRAY_BUFFER);

		// Vertex positions
		bunny.V.resize(3, 0);
		bunny.V_vbo.update(bunny.V);

		// flat normal
		bunny.N_vbo.init(GL_FLOAT, GL_ARRAY_BUFFER);
		bunny.N.resize(3, 0);
		bunny.N_vbo.update(bunny.N);

		// vertex normal
		bunny.NV_vbo.init(GL_FLOAT, GL_ARRAY_BUFFER);
		bunny.NV.resize(3, 0);
		bunny.NV_vbo.update(bunny.NV);

		// Create a new VAO for the bunny. and bind it
		bunny.vao.init();
		bunny.vao.bind();

		// Bind the element buffer, this information will be stored in the current VAO

		// The vertex shader wants the position of the vertices as an input.
		// The following line connects the VBO we defined above with the position "slot"
		// in the vertex shader
		program.bindVertexAttribArray("position", bunny.V_vbo);
		program.bindVertexAttribArray("N", bunny.N_vbo);
		program.bindVertexAttribArray("NV", bunny.NV_vbo);

		// Unbind the VAO when I am done
		bunny.vao.unbind();
	}

	// For the first exercises, 'view' and 'proj' will be the identity matrices
	// However, the 'model' matrix must change for each model in the scene
	Eigen::Matrix4f I = Eigen::Matrix4f::Identity();
	program.bind();
	glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, I.data());
	glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE, I.data());
	glUniform1i(program.uniform("mode"), bunny.mode); // bind mode
	glUniform3f(program.uniform("light_pos"), 1, 1, 1); // add light source
	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback);

	// Register the mouse callback
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window)) {
		// Set the size of the viewport (canvas) to the size of the application window (framebuffer)
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);

		// Bind your program
		program.bind();

		{
			// Bind the VAO for the bunny
			bunny.vao.bind();

			// Model matrix for the bunny

			// Set the uniform value depending on the time difference
			glUniform3f(program.uniform("triangleColor"), 0.0f, 0.0f, 0.0f);

			// Draw the triangles
			switch (bunny.mode) {
			case WIREFRAME: {
				glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				for (int i = 0; i < obj_num; ++i) {
					glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny.mtx_list[i]->data());
					glUniform1i(program.uniform("mode"), bunny.mode); // bind mode
					for (int j = 0; j < bunny.obj_v_num[i]; j += 3) {
						glDrawArrays(GL_LINE_LOOP, bunny.off_num[i] + j, 3);
					}
				}
				break;
			}

			case FLAT: {
				glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				for (int i = 0; i < obj_num; ++i) {
					if (bunny.obj_select[i] == 0) {
						glUniform3f(program.uniform("triangleColor"), 1.0f, 1.0f, 1.0f);
					}
					else {
						glUniform3f(program.uniform("triangleColor"), 1.0f, 0.0f, 0.0f); // red
					}
					glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny.mtx_list[i]->data());
					glUniform1i(program.uniform("mode"), bunny.mode); // bind mode
					glDrawArrays(GL_TRIANGLES, bunny.off_num[i], bunny.obj_v_num[i]);
					glUniform1i(program.uniform("mode"), 0); // bind mode
				}
				break;
			}

			case PHONG: {
				glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				for (int i = 0; i < obj_num; ++i) {
					if (bunny.obj_select[i] == 0) {
						glUniform3f(program.uniform("triangleColor"), 1.0f, 1.0f, 1.0f);
					}
					else {
						glUniform3f(program.uniform("triangleColor"), 1.0f, 0.0f, 0.0f); // red
					}
					glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny.mtx_list[i]->data());
					glUniform1i(program.uniform("mode"), bunny.mode); // bind mode
					glDrawArrays(GL_TRIANGLES, bunny.off_num[i], bunny.obj_v_num[i]);
					glUniform1i(program.uniform("mode"), 0); // bind mode
				}
				break;
			}
			}
		}

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	// Deallocate opengl memory
	program.free();
	bunny.vao.free();
	bunny.V_vbo.free();
	bunny.N_vbo.free();
	bunny.NV_vbo.free();

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
