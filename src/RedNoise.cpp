#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
// Lab 2 Task 4
#include <glm/glm.hpp>
#include "glm/ext.hpp"
// Lab 3 Task 2
#include <CanvasPoint.h>
#include <Colour.h>
// Lab 3 Task 3
#include <CanvasTriangle.h>

#define WIDTH 320
#define HEIGHT 240

// Lab 3 helper function
// Task 2 
std::vector<glm::vec2> interpolateTwoElementValues(glm::vec2 from, glm::vec2 to, float numberOfValues){
	std::vector<glm::vec2> result;
	glm::vec2 stepSize = (to - from)/(numberOfValues-1);
	for (float i = 0; i < numberOfValues; i++){
		result.push_back(from + i*stepSize);
	}
	return result;
}

// Task 4 
CanvasTriangle sortTriangle (CanvasTriangle triangle){
	while (!((triangle.v0().y>triangle.v1().y)&&(triangle.v1().y>triangle.v2().y))){
		if (triangle.v0().y<triangle.v1().y) std::swap(triangle.v0().y,triangle.v1().y);
		if (triangle.v1().y<triangle.v2().y) std::swap(triangle.v1().y,triangle.v2().y);
	}
	CanvasTriangle sortedTriangle = {triangle.v0(), triangle.v1(), triangle.v2()};
	return sortedTriangle;
}

CanvasPoint findMiddlePoint (CanvasTriangle triangle){
	float x0 = triangle.v0().x;
	float x2 = triangle.v2().x;
	float y0 = triangle.v0().y;
	float y1 = triangle.v1().y;
	float y2 = triangle.v2().y;
	float x = x2+(y1-y2)*(x0-x2)/(y0-y2);
	return {x,y1};
}

CanvasTriangle generateRandomTriangle(DrawingWindow &window){
	CanvasPoint v0 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasPoint v1 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasPoint v2 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasTriangle aLovelyTriangle = CanvasTriangle(v0, v1, v2);
	return aLovelyTriangle;
}


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

// Lab 3 Task 2
void drawLine (DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour color){
	std::vector<glm::vec2> interpolated = 
	interpolateTwoElementValues({from.x, from.y},
								{to.x, to.y}, 
								std::max(abs(to.x-from.x), abs(to.y-from.y)));
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + color.blue;
	for (float i = 0.0; i < interpolated.size(); i++) {
		window.setPixelColour(interpolated[i][0], interpolated[i][1], colour);
	}
}

// Lab 3 Task 3
void drawStrokedTriangle (DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color){
	drawLine(window, v0, v1, color);
	drawLine(window, v1, v2, color);
	drawLine(window, v0, v2, color);
}

void drawHalfTriangle (DrawingWindow &window, CanvasPoint h0, CanvasPoint h1, CanvasPoint h2, Colour color){
	std::vector<glm::vec2> h0h1  = interpolateTwoElementValues({h0.x, h0.y}, {h1.x, h1.y}, h1.y-h0.y);
	std::vector<glm::vec2> h0h2  = interpolateTwoElementValues({h0.x, h0.y}, {h2.x, h2.y}, h1.y-h0.y);
	for (float i = 0.0; i < h0h1.size(); i++){
		drawLine(window, {h0h1[i][0], h0h1[i][1]}, {h0h2[i][0], h0h2[i][1]}, color);
	}
}

// Lab 3 Task 4
void drawFilledTriangle (DrawingWindow &window, CanvasTriangle triangle, Colour color){
	Colour white = {255, 255, 255};
	triangle = sortTriangle(triangle);
	CanvasPoint m = findMiddlePoint(triangle);
	drawStrokedTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), white);
	// drawHalfTriangle(window, triangle.v0(), m, triangle.v1(), color);
	drawHalfTriangle(window, triangle.v2(), triangle.v1(), m, color);
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	Colour color = {rand()%255, rand()%255, rand()%255};
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		// Lab 3 Task 3
		else if (event.key.keysym.sym == 'u') {
			Colour color = {rand()%255, rand()%255, rand()%255};
			CanvasTriangle triangle = generateRandomTriangle(window);
			drawStrokedTriangle(window, triangle.v0(),triangle.v1(),triangle.v2(), color);
		}
		else if (event.key.keysym.sym == 'f') {
			CanvasTriangle triangle = generateRandomTriangle(window);
			drawFilledTriangle(window, triangle, color);
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	// Lab 3 Task 2
	// Colour color = {255,255,255};
	// float w = (float)window.width-1;
	// float h = (float)window.height-1;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// Lab 3 Task 2 test
		// drawLine(window, {0.0,0.0}, {w/2,h/2}, color); // A line from the top-left corner of the window to the centre of the window.
		// drawLine(window, {w/2,0.0}, {w/2,h}, color); // A vertical line all the way down the middle of the screen.
		// drawLine(window, {w/3,h/2}, {2*w/3,h/2}, color); // A horizontal line a third the width of the screen, centred both horizontally and vertically.
		// drawLine(window, {w,0.0}, {0.0,h}, color); // A line from the top-right corner to the bottom-left corner.
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
