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
 

void cube_init(Eigen::MatrixXf &cube_V) {
	// init Vertex
	cube_V.resize(36, 3);
	cube_V << -0.5, -0.5, -0.5,
			  -0.5, -0.5,  0.5,
			   0.5, -0.5,  0.5,

			  -0.5, -0.5, -0.5,
			   0.5, -0.5,  0.5,
			   0.5, -0.5, -0.5,

			  -0.5, -0.5, -0.5,
			  -0.5,  0.5, -0.5,
			   0.5,  0.5, -0.5,

			  -0.5, -0.5, -0.5,
			   0.5,  0.5, -0.5,
			   0.5, -0.5, -0.5,

			  -0.5, -0.5, -0.5,
			  -0.5,  0.5, -0.5,
			  -0.5,  0.5,  0.5,

			  -0.5, -0.5, -0.5,
			  -0.5,  0.5,  0.5,
			  -0.5, -0.5,  0.5,

			   0.5,  0.5,  0.5,
			   0.5, -0.5, -0.5,
			   0.5, -0.5,  0.5,

			   0.5,  0.5,  0.5,
			   0.5,  0.5, -0.5,
			   0.5, -0.5, -0.5,

			   0.5,  0.5,  0.5,
			   0.5, -0.5,  0.5,
			  -0.5, -0.5,  0.5,

			   0.5,  0.5,  0.5,
			  -0.5, -0.5,  0.5,
			  -0.5,  0.5,  0.5,

			   0.5,  0.5,  0.5,
			  -0.5,  0.5,  0.5,
			  -0.5,  0.5, -0.5,

			   0.5,  0.5,  0.5,
			  -0.5,  0.5, -0.5,
			   0.5,  0.5, -0.5;
	cube_V.transposeInPlace();
}


////////////////////////////////////////////////////////////////////////////////


// Mesh object, with both CPU data (Eigen::Matrix) and GPU data (the VBOs)
struct Mesh {
	Eigen::MatrixXf V; // mesh vertices [3 x n]

	Eigen::MatrixXf T; // mesh trasformation matrix [4 x j]
	Eigen::VectorXi TI; // mesh transformation index
	std::vector<Eigen::Matrix4f *> mtx_list; // track transform mtx
	std::vector<int> obj_v_num;  // track how many triangle per obj 
	std::vector<int> off_num = {0};  // track offset

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
void load_off(const std::string &filename, Eigen::MatrixXf &V) {
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
	for (int i = 0; i < nf; ++i) {
		V.col(3 * i) = tmpV.col(tmpF(0, i));
		V.col(3 * i + 1) = tmpV.col(tmpF(1, i));
		V.col(3 * i + 2) = tmpV.col(tmpF(2, i));
	}
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
	double highdpi = (double) width / (double) width_window;
	xpos *= highdpi;
	ypos *= highdpi;

	// Convert screen position to the canonical viewing volume
	double xcanonical = ((xpos/double(width))*2)-1;
	double ycanonical = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw

	// TODO: Ray-casting for object selection (Ex.3)
}

void add_off(std::string dir) {
	Eigen::MatrixXf Tmp_V;
	std::string full_dir = DATA_DIR + dir;
	load_off(full_dir , Tmp_V);
	// add offset
	long orig_V_cols = bunny.V.cols();
	bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + Tmp_V.cols()); // resize for cube
	for (int j = 0; j < Tmp_V.cols(); ++j) {
		bunny.V.col(orig_V_cols + j) = Tmp_V.col(j);
	} 

	bunny.V_vbo.update(bunny.V);
	float scale = 0;
	float x_off = - (Tmp_V.row(0).minCoeff() + Tmp_V.row(0).maxCoeff()) / 2;
	float y_off = - (Tmp_V.row(1).minCoeff() + Tmp_V.row(1).maxCoeff()) / 2;
	float z_off = - (Tmp_V.row(2).minCoeff() + Tmp_V.row(2).maxCoeff()) / 2;
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
	scale_mtx(3,3) = 1;
	*trans = scale_mtx * (*trans);
	bunny.mtx_list.push_back(trans);
	bunny.obj_v_num.push_back(Tmp_V.cols());
	bunny.off_num.push_back(bunny.off_num.back() + Tmp_V.cols());
	obj_num++; // increment object count
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS) {
	switch (key) {
			case GLFW_KEY_1: {// add cube
				Eigen::MatrixXf cube_V;
				cube_init(cube_V);
				long orig_V_cols = bunny.V.cols();
				bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + cube_V.cols()); // resize for cube
				for (int j = 0; j < cube_V.cols(); ++j) {
					bunny.V.col(orig_V_cols + j) = cube_V.col(j);
				} 
				bunny.V_vbo.update(bunny.V);
				Eigen::Matrix4f* trans = new Eigen::Matrix4f; // transformation matrix
				*trans = Eigen::Matrix4f::Identity(); // add identity
				bunny.mtx_list.push_back(trans);
				bunny.obj_v_num.push_back(cube_V.cols());
				bunny.off_num.push_back(bunny.off_num.back() + cube_V.cols());
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
				bunny.off_num = {0};
				obj_num = 0; // reset 
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
	GLFWwindow * window = glfwCreateWindow(640, 640, "[Float] Hello World", NULL, NULL);
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

		void main() {
			gl_Position = proj * view * model * vec4(position, 1.0);
		}
	)";

	const GLchar* fragment_shader = R"(
		#version 150 core

		uniform vec3 triangleColor;
		out vec4 outColor;

		void main() {
			outColor = vec4(triangleColor, 1.0);
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
		// Create a new VAO for the bunny. and bind it
		bunny.vao.init();
		bunny.vao.bind();

		// Bind the element buffer, this information will be stored in the current VAO

		// The vertex shader wants the position of the vertices as an input.
		// The following line connects the VBO we defined above with the position "slot"
		// in the vertex shader
		program.bindVertexAttribArray("position", bunny.V_vbo);

		// Unbind the VAO when I am done
		bunny.vao.unbind();
	}

	// For the first exercises, 'view' and 'proj' will be the identity matrices
	// However, the 'model' matrix must change for each model in the scene
	Eigen::Matrix4f I = Eigen::Matrix4f::Identity();
	program.bind();
	glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, I.data());
	glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE, I.data());

	// Save the current time --- it will be used to dynamically change the triangle color
	auto t_start = std::chrono::high_resolution_clock::now();

	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback);

	// Register the mouse callback
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window)) {
		// Set the size of the viewport (canvas) to the size of the application window (framebuffer)
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);

		// Clear the framebuffer
		glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Bind your program
		program.bind();

		{
			// Bind the VAO for the bunny
			bunny.vao.bind();

			// Model matrix for the bunny

			// Set the uniform value depending on the time difference
			auto t_now = std::chrono::high_resolution_clock::now();
			float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
			glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);

			// Draw the triangles
			for (int i = 0; i < obj_num; ++i) {
				glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny.mtx_list[i]->data());
				for (int j = 0; j < bunny.obj_v_num[i]; j+=3) {
					glDrawArrays(GL_LINE_LOOP, bunny.off_num[i] + j, 3);
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

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
