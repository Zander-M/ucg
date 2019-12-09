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
std::vector<Eigen::Matrix4f *> mtx_list; // track transform mtx
std::vector<int> obj_tri_num;  // track how many triangle per obj 
std::vector<int> off_num = {0};  // track offset 

void cube_init(Eigen::MatrixXf &cube_V, Eigen::MatrixXi &cube_F, int len) {
	// init Vertex
	cube_V.resize(8, 3);
	cube_F.resize(12, 3);
	// vertex coord
	cube_V << -0.5, -0.5, -0.5,
			    0.5, -0.5, -0.5,
			    0.5,  0.5, -0.5,
			   -0.5,  0.5, -0.5,
			   -0.5, -0.5,  0.5,
			    0.5, -0.5,  0.5,
			    0.5,  0.5,  0.5,
			   -0.5,  0.5,  0.5;

	// vertex index
	cube_F <<  0 + len, 1 + len, 2 + len,
			   0 + len, 3 + len, 2 + len,
			   0 + len, 1 + len, 5 + len, 
			   0 + len, 4 + len, 5 + len,
			   0 + len, 3 + len, 7 + len,
			   0 + len, 4 + len, 7 + len,
			   6 + len, 5 + len, 1 + len,
			   6 + len, 2 + len, 1 + len,
			   6 + len, 5 + len, 4 + len,
			   6 + len, 7 + len, 4 + len,
			   6 + len, 2 + len, 3 + len,
			   6 + len, 7 + len, 3 + len;
	cube_V.transposeInPlace();
	cube_F.transposeInPlace();
}


////////////////////////////////////////////////////////////////////////////////


// Mesh object, with both CPU data (Eigen::Matrix) and GPU data (the VBOs)
struct Mesh {
	Eigen::MatrixXf V; // mesh vertices [3 x n]
	Eigen::MatrixXi F; // mesh triangles [3 x m]
	Eigen::MatrixXf T; // mesh trasformation matrix [4 x j]
	Eigen::VectorXi TI; // mesh transformation index

	// VBO storing vertex position attributes
	VertexBufferObject V_vbo;

	// VBO storing vertex indices (element buffer)
	VertexBufferObject F_vbo;

	// VCO storing transformations

	// VCO storing transformation index buffer

	// VAO storing the layout of the shader program for the object 'bunny'
	VertexArrayObject vao;
};

Mesh bunny;

////////////////////////////////////////////////////////////////////////////////

// Read a triangle mesh from an off file
void load_off(const std::string &filename, Eigen::MatrixXf &V, Eigen::MatrixXi &F) {
	std::ifstream in(filename);
	std::string token;
	in >> token;
	int nv, nf, ne;
	in >> nv >> nf >> ne;
	V.resize(3, nv);
	F.resize(3, nf);
	for (int i = 0; i < nv; ++i) {
		in >> V(0, i) >> V(1, i) >> V(2, i);
	}
	for (int i = 0; i < nf; ++i) {
		int s;
		in >> s >> F(0, i) >> F(1, i) >> F(2, i);
		assert(s == 3);
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
	Eigen::MatrixXi Tmp_F;
	std::string full_dir = DATA_DIR + dir;
	load_off(full_dir , Tmp_V, Tmp_F);
	// add offset
	Eigen::MatrixXi Offset;
	Offset = Eigen::MatrixXi::Constant(Tmp_F.rows(), Tmp_F.cols(), bunny.V.cols());
	Tmp_F = Tmp_F + Offset; 
	long orig_V_cols = bunny.V.cols();
	long orig_F_cols = bunny.F.cols();
	bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + Tmp_V.cols()); // resize for cube
	for (int j = 0; j < Tmp_V.cols(); ++j) {
		bunny.V.col(orig_V_cols + j) = Tmp_V.col(j);
	} 
	bunny.F.conservativeResize(bunny.F.rows(), bunny.F.cols() + Tmp_F.cols()); // resize for cube
	for (int j = 0; j < Tmp_F.cols(); ++j) {
		bunny.F.col(orig_F_cols + j) = Tmp_F.col(j);
	} 
	bunny.V_vbo.update(bunny.V);
	bunny.F_vbo.update(bunny.F);
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
	mtx_list.push_back(trans);
	obj_tri_num.push_back(Tmp_F.cols());
	off_num.push_back(off_num.back() + 3 * Tmp_F.cols());
	obj_num++; // increment object count
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS) {
	switch (key) {
			case GLFW_KEY_1: {// add cube
				Eigen::MatrixXf cube_V;
				Eigen::MatrixXi cube_F;
				cube_init(cube_V, cube_F, bunny.V.cols());
				long orig_V_cols = bunny.V.cols();
				long orig_F_cols = bunny.F.cols();
				bunny.V.conservativeResize(bunny.V.rows(), bunny.V.cols() + 8); // resize for cube
				for (int j = 0; j < 8; ++j) {
					bunny.V.col(orig_V_cols + j) = cube_V.col(j);
				} 
				bunny.F.conservativeResize(bunny.F.rows(), bunny.F.cols() + 12); // resize for cube
				for (int j = 0; j < 12; ++j) {
					bunny.F.col(orig_F_cols + j) = cube_F.col(j);
				} 
				bunny.V_vbo.update(bunny.V);
				bunny.F_vbo.update(bunny.F);
				Eigen::Matrix4f* trans = new Eigen::Matrix4f; // transformation matrix
				*trans = Eigen::Matrix4f::Identity(); // add identity
				mtx_list.push_back(trans);
				obj_tri_num.push_back(cube_F.cols());
				off_num.push_back(off_num.back() + 3 * cube_F.cols());
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
				bunny.F.resize(3, 0);
				bunny.V_vbo.update(bunny.V);
				bunny.F_vbo.update(bunny.F);
				for (int i = 0; i < obj_num; ++i) {
					delete mtx_list[i];
				}
				mtx_list.clear();
				obj_tri_num.clear();
				off_num = {0};
				obj_num = 0; // reset 
				break;
			}

			case GLFW_KEY_P: {
				printf("rows: %ld\n", bunny.V.cols());
				printf("obj_count: %d\n", obj_num);
				std::cout << "Offset:" << std::endl;
				for (int i = 0; i < off_num.size(); ++i)
					std::cout << off_num[i] << " ";
				std::cout << std::endl;
				std::cout << "Tri num:" << std::endl;
				for (int i = 0; i < obj_tri_num.size(); ++i)
					std::cout << obj_tri_num[i] << " ";
				std::cout << std::endl;
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
		bunny.F_vbo.init(GL_UNSIGNED_INT, GL_ELEMENT_ARRAY_BUFFER);

		// Vertex positions
		bunny.V.resize(3, 0);
		// bunny.V <<
		// 	0, 0.5, -0.5,
		// 	0.5, -0.5, -0.5,
		// 	0, 0, 0;
		bunny.V_vbo.update(bunny.V);

		// Triangle indices
		bunny.F.resize(3, 0);
		// bunny.F << 0, 1, 2;
		bunny.F_vbo.update(bunny.F);

		// Create a new VAO for the bunny. and bind it
		bunny.vao.init();
		bunny.vao.bind();

		// Bind the element buffer, this information will be stored in the current VAO
		bunny.F_vbo.bind();

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
				glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, mtx_list[i]->data());
				// for (int j = 0; j < obj_tri_num[i]; ++j) {
				// 	glDrawArrays(GL_TRIANGLES, off_num[i] + 3*j, 3);
				// }
				// for (int j = 0; j < obj_tri_num[i]; ++j) {
				// 	glDrawArrays(GL_TRIANGLES, off_num[i] + 3 * j, 3);
				// }
				glDrawElements(GL_LINE_LOOP, off_num[i+1], bunny.F_vbo.scalar_type, NULL);
			}
			// for (int i = 0; i < bunny.F.cols()/3; ++i) {
				// glDrawArrays(GL_LINE_LOOP, i * 3, 3);
			// }
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
	bunny.F_vbo.free();

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
