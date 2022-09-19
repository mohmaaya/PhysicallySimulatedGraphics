

#include <thread>



//#include <include/GL/glut.h>


#include <glad/include/glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Gui.h"
#include "Simulator.h"
#include "FluidSim.h"

#include "meshFire.h"
#include "mesh22.h"
#include "shader.h"
#include "model.h"
#include "camera.h"
#include "animator.h"

#include <math.h>

#include "fire_and_air.h"

#include "gridHelper.h"

#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>




	void framebuffer_size_callback(GLFWwindow* window, int width, int height);
	void mouse_callback(GLFWwindow* window, double xpos, double ypos);
	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
	void processInput(GLFWwindow* window);

	// settings
	const unsigned int SCR_WIDTH = 800;
	const unsigned int SCR_HEIGHT = 600;

	// camera
	Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
	float lastX = SCR_WIDTH / 2.0f;
	float lastY = SCR_HEIGHT / 2.0f;
	bool firstMouse = true;

	// timing
	float deltaTime = 0.f;
	float lastFrame = 0.f;

	glm::vec3 offset = glm::vec3(0, 0, 0);



	int main(int argc, char* argv[]) {


			// glfw: initialize and configure
		// ------------------------------
		glfwInit();
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif

		// glfw window creation
		// --------------------
		GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
		if (window == NULL)
		{
			std::cout << "Failed to create GLFW window" << std::endl;
			glfwTerminate();
			return -1;
		}
		glfwMakeContextCurrent(window);
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
		glfwSetCursorPosCallback(window, mouse_callback);
		glfwSetScrollCallback(window, scroll_callback);

		// tell GLFW to capture our mouse
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

		// glad: load all OpenGL function pointers
		// ---------------------------------------
		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
			std::cout << "Failed to initialize GLAD" << std::endl;
			return -1;
		}

		// tell stb_image.h to flip loaded texture's on the y-axis (before loading model).
		stbi_set_flip_vertically_on_load(true);

		// configure global opengl state
		// -----------------------------
		glEnable(GL_DEPTH_TEST);

		// build and compile shaders
		// -------------------------
		Shader ourShader("C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/basic_vs.glsl", "C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/basic_fs.glsl");
		Shader fireairShader("C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/fire_air_vs.glsl", "C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/fire_air_fs.glsl");


		// load models
		// -----------


	

		Model ourModel("C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/avatar1.glb");
		Animation danceAnimation("C:/Users/moh/Documents/FAU/PBS/Projekt/version01/pbs-ws21/6_avatar/avatar1.glb", &ourModel);
		Animator animator(&danceAnimation);

		float sizegrid = 32 * 32 * 32 * sizeof(float);

		Fire_and_Air fireandair(1);
		GridData gridfire;

		gridfire.initMesh(glm::vec3(-0.5, 1., 0.2), 0.05f, glm::vec3(31, 31, 31));
		//gridfire.initMesh(glm::vec3(-0.4, 1.1, 0.2), 0.009f, glm::vec3(31, 31, 31));

		MeshFire fireMesh(gridfire.vertices, gridfire.indices);


// draw in wireframe
//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

// render loop
// -----------

//std::vector<float> test;
//test.assign(32 * 32 * 32, 0.5f);
//
//for (int i = 0; i < 10; i++)
//	fireandair.fire->T[i] = i;


		unsigned int fireBuffer;
		glGenBuffers(1, &fireBuffer);
		glBindBuffer(GL_ARRAY_BUFFER, fireBuffer);

		if (fireandair.showFire == 1)
		glBufferData(GL_ARRAY_BUFFER, sizegrid, fireandair.fire->T, GL_STATIC_DRAW);

		if (fireandair.showAir == 1)
			glBufferData(GL_ARRAY_BUFFER, sizegrid, fireandair.air->d, GL_STATIC_DRAW);

		unsigned int fireAttrib = glGetAttribLocation(fireairShader.ID, "fire");

		glEnableVertexAttribArray(fireAttrib);
		glVertexAttribPointer(fireAttrib, 1, GL_FLOAT, GL_FALSE, 0, 0);

		


while (!glfwWindowShouldClose(window))
{
	// per-frame time logic
	// --------------------
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;

	// input
	// -----
	processInput(window);


	

	if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
		fireandair.key_func('x', 1, 1);
	if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
		fireandair.key_func('y', 1, 1);
	if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
		fireandair.key_func('n', 1, 1);
	if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
		fireandair.key_func('z', 1, 1);
	if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
		fireandair.key_func('f', 1, 1);
	if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS)
		fireandair.key_func('i', 1, 1);



	animator.UpdateAnimation(deltaTime);

	// render
	// ------
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// don't forget to enable shader before setting uniforms
	ourShader.use();

	// view/projection transformations
	glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
	glm::mat4 view = camera.GetViewMatrix();
	ourShader.setMat4("projection", projection);
	ourShader.setMat4("view", view);


	auto transforms = animator.GetFinalBoneMatrices();
	for (int i = 0; i < transforms.size(); ++i)
		ourShader.setMat4("finalBonesMatrices[" + std::to_string(i) + "]", transforms[i]);

	// render the loaded model
	glm::mat4 model = glm::mat4(1.0f);
	model = glm::translate(model, glm::vec3(0.0f, -0.4f, 0.0f)); // translate it down so it's at the center of the scene
	model = glm::scale(model, glm::vec3(.05f, .05f, .05f));	// it's a bit too big for our scene, so scale it down
	ourShader.setMat4("model", model);
	ourModel.Draw(ourShader);


	fireairShader.use();


	glBindBuffer(GL_ARRAY_BUFFER, fireBuffer);


	if (fireandair.showFire == 1) 
		glBufferData(GL_ARRAY_BUFFER, sizegrid, fireandair.fire->T, GL_STATIC_DRAW);
		

	if (fireandair.showAir == 1)
		glBufferData(GL_ARRAY_BUFFER, sizegrid, fireandair.air->d, GL_STATIC_DRAW);


	if (fireandair.showFire == 1)
	fireairShader.setVec3("colorFire", glm::vec3(1, 0.2, 0));

	if (fireandair.showAir == 1)
	fireairShader.setVec3("colorFire", glm::vec3(0.65, 0.9, 1));
	fireairShader.setVec4("colorDefault", glm::vec4(1, 1, 1, 0));
	fireairShader.setFloat("threshhold", 0.1f);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glDrawElements(GL_TRIANGLES, testData.triangleCount * 3, GL_UNSIGNED_INT, 0);



			fireandair.sim_fire_and_air();
			
			fireairShader.setMat4("projection", projection);
			fireairShader.setMat4("view", view);
			fireairShader.setVec3("offset", offset);
			
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDisable(GL_DEPTH_TEST);

			if(fireandair.showFire == 1)
			fireMesh.Draw(fireandair.fire->T, sizegrid); //fire.size

			if (fireandair.showAir == 1)
				fireMesh.Draw(fireandair.air->d, sizegrid); //fire.size

			glEnable(GL_DEPTH_TEST);
			glDisable(GL_BLEND);
			
			glDisable(GL_BLEND);

			// glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
			// -------------------------------------------------------------------------------
			glfwSwapBuffers(window);
			glfwPollEvents();
		}

		// glfw: terminate, clearing all previously allocated GLFW resources.
		// ------------------------------------------------------------------
		glfwTerminate();

		
		return 0;



	}




	// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
	// ---------------------------------------------------------------------------------------------------------
	void processInput(GLFWwindow* window)
	{
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(window, true);

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
			camera.ProcessKeyboard(FORWARD, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
			camera.ProcessKeyboard(BACKWARD, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
			camera.ProcessKeyboard(LEFT, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
			camera.ProcessKeyboard(RIGHT, deltaTime);
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
			offset.x += 0.01;
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
			offset.x -= 0.01;
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
			offset.y+= 0.1;
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
			offset.y -= 0.1;
		if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
			offset.z += 0.1;
		if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
			offset.z -= 0.1;

	}

	// glfw: whenever the window size changed (by OS or user resize) this callback function executes
	// ---------------------------------------------------------------------------------------------
	void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		// make sure the viewport matches the new window dimensions; note that width and 
		// height will be significantly larger than specified on retina displays.
		glViewport(0, 0, width, height);
	}

	// glfw: whenever the mouse moves, this callback is called
	// -------------------------------------------------------
	void mouse_callback(GLFWwindow* window, double xpos, double ypos)
	{
		if (firstMouse)
		{
			lastX = xpos;
			lastY = ypos;
			firstMouse = false;
		}

		float xoffset = xpos - lastX;
		float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

		lastX = xpos;
		lastY = ypos;

		camera.ProcessMouseMovement(xoffset, yoffset);
	}

	// glfw: whenever the mouse scroll wheel scrolls, this callback is called
	// ----------------------------------------------------------------------
	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		camera.ProcessMouseScroll(yoffset);
	}


