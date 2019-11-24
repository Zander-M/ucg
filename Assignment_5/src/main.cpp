// This example is heavily based on the tutorial at https://open.gl

////////////////////////////////////////////////////////////////////////////////
// OpenGL Helpers to reduce the clutter
#include "helpers.h"
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
// Linear Algebra Library
#include <Eigen/Dense>
// Timer
#include <chrono>
////////////////////////////////////////////////////////////////////////////////

enum MODE { // for mode
    NORMAL = 0,
    INSERTION,
    TRANSLATION,
    DELETE,
    COLOR
};

enum V_COLOR { // for vertex color
    COLOR_1 = 0,
    COLOR_2,
    COLOR_3,
    COLOR_4,
    COLOR_5,
    COLOR_6,
    COLOR_7,
    COLOR_8,
    COLOR_9
};

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_color; // VBO for color
MODE mode = NORMAL;
V_COLOR color = COLOR_1;

int t_count = 0;     // count number of triangles
int t_selected = -1; // index of triangle selected
int t_prev = -1;
int v_selected = -1; // index of vertex selected
int active_vtx =
    0; // for insertion mode, keep track of how many vertices have been placed
int trans_mode = 0; // indicator of whether in translation mode

// Contains the vertex positions
Eigen::MatrixXf view(4, 4);      // transformation matrix
Eigen::MatrixXf V(2, 3);          // Vertex Matrix
Eigen::MatrixXf Color(3, 3);      // for RGB color
Eigen::MatrixXf Prev_Color(3, 3); // previous color before selection.
Eigen::Vector2f Click_Pos;        // click position
Eigen::Vector2f Prev_Pos;        // previous position
Eigen::Vector2f offset;           // movement offset

//////// HELPER FUNCTIONS ////////

// return 0if a point intersects with a triangle and set t_selected, else
// -1.
int point2triangle(Eigen::Vector2f pt, Eigen::MatrixXf vertex, int t_cnt) {
    for (int i = t_cnt - 1; i >= 0; --i) { // index start from 0
        Eigen::Vector2f a = V.col(3 * i + 2) - V.col(3 * i);
        Eigen::Vector2f b = V.col(3 * i + 1) - V.col(3 * i);
        Eigen::Matrix2f mtx;
        mtx.col(0) = a;
        mtx.col(1) = b;
        Eigen::Vector2f sol =
            mtx.colPivHouseholderQr().solve(pt - V.col(3 * i));
        if (sol(0) > 0 && sol(1) > 0 && sol(0) + sol(1) < 1) {
            t_selected = i;
            return 0;
        }
    }
    t_selected = -1;
    return -1;
}

// return index of the closest vertex, if none return -1
int point2vertex(Eigen::Vector2f pt, Eigen::MatrixXf vertex) {
    double dis = 3; // set click distance
    int ind = -1;
    Eigen::Vector2f diff;
    for (int i = 0; i < V.cols(); ++i) {
        diff = pt - V.col(i);
        if (diff.norm() < dis && diff.norm() < 0.2) {
            ind = i;
            dis = diff.norm();
        }
    }
    v_selected = ind;
    printf("vertex: %d\n", ind);
    return 0;
}

/////// CALLBACKS ////////
void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
    // Get viewport size (canvas in number of pixels)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    // Get the size of the window (may be different than the canvas size on
    // retina displays)
    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Deduce position of the mouse in the viewport
    double highdpi = (double)width / (double)width_window;
    xpos *= highdpi;
    ypos *= highdpi;

    // Convert screen position to world coordinates
	Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
    Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    Eigen::Vector4f p_world = view.inverse()*p_canonical;
    Click_Pos << p_world(0), p_world(1); // update click pos

    // Update the position of the first vertex if the left button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        switch (mode) {
        case NORMAL:
            break;
        case INSERTION: {
            printf("INSERTION\n");
            active_vtx++;          // increment vertex
            if (active_vtx == 1) { // if it's ths first vertex
                V.conservativeResize(V.rows(),
                                     V.cols() + 3); // add three new nodes
                Color.conservativeResize(
                    Color.rows(), Color.cols() + 3); // add three new nodes
                V.col(V.cols() - 3) << Click_Pos;
                V.col(V.cols() - 2) << Click_Pos;
                V.col(V.cols() - 1) << Click_Pos;
                Color.col(Color.cols() - 3) << 0, 0, 0; // initialize as black
                Color.col(Color.cols() - 2) << 0, 0, 0;
                Color.col(Color.cols() - 1) << 0, 0, 0;
            } else if (active_vtx == 2) {
                V.col(V.cols() - 2) << Click_Pos;
                V.col(V.cols() - 1) << Click_Pos;
            } else if (active_vtx == 3) {
                V.col(V.cols() - 1) << Click_Pos;
                Color.col(Color.cols() - 3) << 1, 0, 0; // set color to red
                Color.col(Color.cols() - 2) << 1, 0, 0;
                Color.col(Color.cols() - 1) << 1, 0, 0;
                active_vtx = 0; // reset
                t_count++;      // add one
            }
            break;
        }
        case TRANSLATION: {
            printf("TRANSLATION\n");
            point2triangle(Click_Pos, V, t_count);
            if (t_selected != -1) {
                Prev_Color.col(0) << Color.col(3 * t_selected);
                Prev_Color.col(1) << Color.col(3 * t_selected + 1);
                Prev_Color.col(2) << Color.col(3 * t_selected + 2);
                Color.col(3 * t_selected) << 0, 0,
                    1; // change color to blue
                Color.col(3 * t_selected + 1) << 0, 0, 1;
                Color.col(3 * t_selected + 2) << 0, 0, 1;
                Prev_Pos << p_world(0), p_world(1);
                trans_mode = 1;
            }
            break;
        }
        case DELETE: {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                point2triangle(Click_Pos, V, t_count);
                if (t_selected != -1) {
                    for (int i = t_selected * 3; i < V.cols() - 3; ++i) {
                        V.col(i) << V.col(i + 3);
                        Color.col(i) << Color.col(i + 3);
                    }
                    V.conservativeResize(V.rows(), V.cols() - 3);
                    Color.conservativeResize(Color.rows(), Color.cols() - 3);
                    t_count--;
                    t_selected = -1;
                }
            }
            break;
        }
        case COLOR: {
			printf("COLOR\n");
            point2vertex(Click_Pos, V);
            break;
        }
        default: {
            printf("Shouldn't get here!\n");
            break;
        }
        }
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        if (mode == TRANSLATION && trans_mode == 1) {
            // offset << xworld - Click_Pos(0), yworld - Click_Pos(1);
            // V.col(3 * t_selected ) += offset;
            // V.col(3 * t_selected + 1) += offset;
            // V.col(3 * t_selected + 2) += offset;
            // restore color
            Color.col(3 * t_selected) << Prev_Color.col(0);
            Color.col(3 * t_selected + 1) << Prev_Color.col(1);
            Color.col(3 * t_selected + 2) << Prev_Color.col(2);
            t_prev = t_selected;
            t_selected = -1; // reset
            trans_mode = 0;
        }
    }
    // Upload the change to the GPU
    VBO.update(V);
    VBO_color.update(Color);
}

void cursor_move_callback(GLFWwindow *window, double xpos, double ypos) {
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    // D	educe position of the mouse in the viewport
    double highdpi = (double)width / (double)width_window;
    xpos *= highdpi;
    ypos *= highdpi;

	Eigen::Vector4f p_screen(xpos,height-1-ypos,0,1);
    Eigen::Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    Eigen::Vector4f p_world = view.inverse()*p_canonical;
    // NOTE: y axis is flipped in glfw
    // printf("%f %f\n", p_world(0), p_world(1)); // test
    switch (mode) {
    case NORMAL:
        break;
    case INSERTION: {
        if (active_vtx == 1) {
            V.col(V.cols() - 2) << p_world(0), p_world(1);
            V.col(V.cols() - 1) << p_world(0), p_world(1);
        } else if (active_vtx == 2) {
            V.col(V.cols() - 1) << p_world(0), p_world(1);
        }
        break;
    }
    case TRANSLATION: {
        if (t_selected != -1 && trans_mode == 1) {
            offset << p_world(0) - Prev_Pos(0), p_world(1) - Prev_Pos(1);
            V.col(t_selected * 3) += offset;
            V.col(t_selected * 3 + 1) += offset;
            V.col(t_selected * 3 + 2) += offset;
            Prev_Pos << p_world(0), p_world(1);
        }
        break;
    }
    default:
        break;
    }
    VBO.update(V);
    VBO_color.update(Color);
}

void key_callback(GLFWwindow *window, int key, int scancode, int action,
                  int mods) {
    // Update the position of the first vertex if the keys 1,2, or 3 are
    // pressed
    if (action == GLFW_PRESS) {
    switch (key) {
    case GLFW_KEY_1: // color 1
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 1, 0, 0; // red
        }
        break;
    case GLFW_KEY_2: // color 2
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 0, 1, 0; // red
        }
        break;
    case GLFW_KEY_3: // color 3
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 0, 0, 1; // red
        break;
    case GLFW_KEY_4: // color 4
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 1, 1, 0; // red
        }
        break;
    case GLFW_KEY_5: // color 5
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 1, 0, 1; // red
        }
        break;
    case GLFW_KEY_6: // color 6
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 0, 1, 1; // red
        }
        break;
    case GLFW_KEY_7: // color 7
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 1, 1, 1; // red
        }
        break;
    case GLFW_KEY_8: // color 8
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 0.5, 0.5, 0; // red
        }
        break;
    case GLFW_KEY_9: // color 9
        if (mode == COLOR && v_selected != -1) {
            Color.col(v_selected) << 1, 0.5, 0.5; // red
        }
        break;
    case GLFW_KEY_I:
        mode = INSERTION;
        printf("INSERTION\n");
        break;
    case GLFW_KEY_O:
        if (mode == INSERTION && active_vtx != 0) {
            V.conservativeResize(V.rows(),
                                 V.cols() - 3); // if previous mode is insertion
                                                // and triangle is not drawn
            Color.conservativeResize(Color.rows(), Color.cols() - 3);
        }
        mode = TRANSLATION;
        printf("TRANSLATION\n");
        break;
    case GLFW_KEY_P:
        if (mode == INSERTION && active_vtx != 0) {
            V.conservativeResize(V.rows(),
                                 V.cols() - 3); // if previous mode is insertion
                                                // and triangle is not drawn
            Color.conservativeResize(Color.rows(), Color.cols() - 3);
        }
        mode = DELETE;
        printf("DELETE\n");
        break;
    case GLFW_KEY_C:
        if (mode == INSERTION && active_vtx != 0) {
            V.conservativeResize(V.rows(),
                                 V.cols() - 3); // if previous mode is insertion
                                                // and triangle is not drawn
            Color.conservativeResize(Color.rows(), Color.cols() - 3);
        }
        mode = COLOR;
        printf("COLOR\n");
        break;
    case GLFW_KEY_Q:
        if (mode == INSERTION && active_vtx != 0) {
            V.conservativeResize(V.rows(),
                                 V.cols() - 3); // if previous mode is insertion
                                                // and triangle is not drawn
            Color.conservativeResize(Color.rows(), Color.cols() - 3);
        }
        mode = NORMAL;
        printf("NORMAL\n");
        break;
    case GLFW_KEY_H:
        break;
    case GLFW_KEY_J:
        break;
    case GLFW_KEY_K:
        break;
    case GLFW_KEY_L:
        break;
    case GLFW_KEY_M:
        printf("%d\n", mode);
        break;
	case GLFW_KEY_F:
		printf("t_cound: %d\n", t_count);
		printf("v_count: %d\n", V.cols());
		printf("t_selected: %d\n", t_selected);
		printf("v_selected: %d\n", v_selected);
        printf("active_v: %d\n", active_vtx);
        printf("t_prev: %d\n", t_prev);
		break;
    default:
        break;
    }
    }
    active_vtx = 0; // set active_back to 0

    // Upload the change to the GPU
    VBO.update(V);
    VBO_color.update(Color);
}
}


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
    GLFWwindow *window =
        glfwCreateWindow(640, 480, "[Float] Hello World", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Load OpenGL and its extensions
    if (!gladLoadGL()) {
        printf("Failed to load OpenGL and its extensions");
        return (-1);
    }
    printf("OpenGL Version %d.%d loaded", GLVersion.major, GLVersion.minor);

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char *)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n",
           (const char *)glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the
    // vertex attributes are stored in a Vertex Buffer Object (or VBO). This
    // means that the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();

    // Initialize the VBO with the vertices data
    // A VBO is a data container that lives in the GPU memory
    VBO.init();
    V.resize(2, 3);
    V << 0, 0.5, -0.5, 0.5, -0.5, -0.5;
    VBO.update(V);

    VBO_color.init();
    Color.resize(3, 3);
    Color << 1, 1, 1, 0, 0, 0, 0, 0, 0; // red
    VBO_color.update(Color);

	t_count = 1; // update t_count

    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar *vertex_shader = R"(
		#version 150 core
		uniform mat4 view;
		in vec2 position;
		in vec3 v_color;
		out vec3 v_colorOut;

		void main() {
			v_colorOut  = v_color;
			gl_Position = view * vec4(position, 0.0, 1.0);
		}
	)";

    const GLchar *fragment_shader = R"(
		#version 150 core

		in vec3 v_colorOut;
		out vec4 outColor;

		void main() {
				outColor = vec4(v_colorOut, 1.0);
		}
	)";

    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called
    // outColor is the one that we want in the fragment buffer (and thus on
    // screen)
    program.init(vertex_shader, fragment_shader, "outColor");
    program.bind();

    // The vertex shader wants the position of the vertices as an input.
    // The following line connects the VBO we defined above with the
    // position "slot" in the vertex shader
    // program.bindVertexAttribArray("trans");
    program.bindVertexAttribArray("position", VBO);
    program.bindVertexAttribArray("v_color", VBO_color);
    glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, view.data());
    // glUniform3f(program.uniform("triangleColor"), 1, 0, 0); // red
    // Save the current time --- it will be used to dynamically change the
    // triangle color
    auto t_start = std::chrono::high_resolution_clock::now();

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    glfwSetCursorPosCallback(window, cursor_move_callback);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        // Set the size of the viewport (canvas) to the size of the
        // application window (framebuffer)
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        float a_ratio = (float)height / (float)width;
        view << 
		a_ratio, 0, 0, 0, // init view matrix, fixed width
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE,
                           view.data()); // bind uniform

        // Bind your VAO (not necessary if you have only one)
        VAO.bind();

        // Bind your program
        program.bind();

        // Set the uniform value depending on the time difference
        // auto t_now = std::chrono::high_resolution_clock::now();
        // float time =
        // std::chrono::duration_cast<std::chrono::duration<float>>(t_now -
        // t_start).count(); glUniform3f(program.uniform("triangleColor"),
        // (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);
        //
        // Clear the framebuffer
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw a triangle
		for (int i = 0; i < t_count; ++i) {
			glDrawArrays(GL_TRIANGLES, i*3, 3);
		}
        if (mode == INSERTION) {
                glDrawArrays(GL_LINE_LOOP, V.cols() - 3, active_vtx + 1);
        }
        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();
    VBO.free();
    VBO_color.free();
    VAO.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
