#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
// Lab 2 Task 4
#include <glm/glm.hpp>
#include "glm/ext.hpp"

#define WIDTH 320
#define HEIGHT 240

// Lab 2 Task 2
std::vector<float> interpolateSingleFloats(float from, float to, float numberOfValues) {
	std::vector<float> result;
	float stepSize = (to - from)/(numberOfValues-1);
	for (int i = 0; i < numberOfValues; i ++) {
		result.push_back(from + i*stepSize);
	}
	return result;
}

// Lab 2 Task 4
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, float numberOfValues){
	std::vector<glm::vec3> result;
	glm::vec3 stepSize = (to - from)/(numberOfValues-1);
	for (float i = 0; i < numberOfValues; i++){ // i must be float. Can't multiply int with float in c++.
		result.push_back(from + i*stepSize);
	}
	return result;
}

// Lab 2 Task 3
// void draw(DrawingWindow &window) {
// 	window.clearPixels();
// 	std::vector<float> colorIndex = interpolateSingleFloats(255, 0, window.width);
// 	for (size_t y = 0; y < window.height; y++) {
// 		for (size_t x = 0; x < window.width; x++) {
// 			float color = colorIndex[x];
// 			uint32_t colour = (255 << 24) + (int(color) << 16) + (int(color) << 8) + int(color);
// 			window.setPixelColour(x, y, colour);
// 		}
// 	}
// }

// Lab 2 Task 5
void draw(DrawingWindow &window) {
	window.clearPixels();
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	std::vector<glm::vec3> y_axisL = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
	std::vector<glm::vec3> y_axisR = interpolateThreeElementValues(topRight, bottomRight, window.height);
	for (size_t j = 0; j < window.height; j ++) {
		std::vector<glm::vec3> x_axis = interpolateThreeElementValues(y_axisL[j], y_axisR[j], window.width);
		for (size_t i = 0; i < window.width; i ++){
			uint32_t colour = (255 << 24) + (int(x_axis[i][0]) << 16) + (int(x_axis[i][1]) << 8) + int(x_axis[i][2]);
			window.setPixelColour(i, j, colour);
		}
	}
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	// Lab 2 Task 2 test
	std::vector<float> result;
	result = interpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	// Lab 2 Task 4 test
	glm::vec3 from(1.0, 4.0, 9.2);
	glm::vec3 to(4.0, 1.0, 9.8);
	std::vector<glm::vec3> result2;
	result2 = interpolateThreeElementValues(from, to, 4);
	for(size_t i=0; i<result2.size(); i++) std::cout << glm::to_string(result2[i]) << " \n ";
	std::cout << std::endl;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
