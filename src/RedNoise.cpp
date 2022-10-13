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
// Lab 3 Task 5
#include <TextureMap.h>

#define WIDTH 320
#define HEIGHT 240

// Lab 3 helper function
// Task 2 
std::vector<glm::vec2> interpolateTwoElementValues(glm::vec2 from, glm::vec2 to, float numberOfValues){
	std::vector<glm::vec2> result;
	glm::vec2 stepSize = (to-from)/(numberOfValues-1);
	for (float i = 0; i < numberOfValues; i++){
		result.push_back(from + i*stepSize);
	}
	return result;
}

// Task 4 
CanvasTriangle sortTriangle (CanvasTriangle triangle){
	while (!((triangle.v0().y<triangle.v1().y)&&(triangle.v1().y<triangle.v2().y))){
		if (triangle.v0().y>triangle.v1().y) std::swap(triangle.v0(),triangle.v1());
		if (triangle.v1().y>triangle.v2().y) std::swap(triangle.v1(), triangle.v2());
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
	return CanvasPoint{x,y1};
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
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;

	for (float i = 0.0; i < numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		uint32_t c = (255 << 24) + (color.red << 16) + (color.green << 8) + color.blue;
		window.setPixelColour(round(x), round(y), c);
	}

}

// Lab 3 Task 3
void drawStrokedTriangle (DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color){
	drawLine(window, v0, v1, color);
	drawLine(window, v1, v2, color);
	drawLine(window, v0, v2, color);
}

// Task 4 helper function
void drawHalfTriangle (DrawingWindow &window, CanvasPoint h0, CanvasPoint h1, CanvasPoint h2, Colour color){
	glm::vec2 from = {h0.x, h0.y};
	glm::vec2 to1 = {h1.x, h1.y};
	glm::vec2 to2 = {h2.x, h2.y};
	float numberOfValues = abs(h1.y - h0.y)+1;
	std::vector<glm::vec2> h0h1  = interpolateTwoElementValues(from, to1, numberOfValues);
	std::vector<glm::vec2> h0h2  = interpolateTwoElementValues(from, to2, numberOfValues);
	for (float i = 0.0; i < numberOfValues; i++){		
		drawLine(window, {round(h0h1[i][0]), h0h1[i][1]}, {h0h2[i][0], h0h2[i][1]}, color);
	}
}

// Lab 3 Task 4
void drawFilledTriangle (DrawingWindow &window, CanvasTriangle triangle, Colour color){
	Colour white = {255, 255, 255};
	triangle = sortTriangle(triangle);
	CanvasPoint m = findMiddlePoint(triangle);
	drawHalfTriangle(window, triangle.v0(), m, triangle.v1(), color);
	drawHalfTriangle(window, triangle.v2(), m, triangle.v1(), color);
	drawStrokedTriangle(window, triangle.v0(), triangle.v1(), triangle.v2(), white);
}

// Task 5 helper function
// TODO: solve the skipped pixels.
void drawHalfTextureTriangle(DrawingWindow &window, CanvasPoint t0, CanvasPoint t1, CanvasPoint t2, TextureMap textureMap, glm::mat3x3 affine){
	float h = abs(t0.y-t1.y);
	std::vector<glm::vec2> t0t1 = interpolateTwoElementValues({t0.x,t0.y},{t1.x,t1.y},h+1); // h+1 to solve skipped lines
	std::vector<glm::vec2> t0t2 = interpolateTwoElementValues({t0.x,t0.y},{t2.x,t2.y},h+1);
	for (float i=0.0; i < h+1; i++){
		std::vector<glm::vec2> line = interpolateTwoElementValues({t0t1[i].x,t0t1[i].y}, {t0t2[i].x,t0t2[i].y}, abs(t0t1[i].x-t0t2[i].x));
		for (float j=0.0; j<line.size(); j++){
			float x = line[j][0];
			float y = line[j][1];
			glm::mat3x3 canvas = glm::mat3{{x,0,0},{y,0,0},{1,0,0}};
			glm::mat3x3 texture = canvas*affine;
			int tx = texture[0][0];
			int ty = texture[1][0];
			uint32_t color = textureMap.pixels[round(ty*textureMap.width+tx)];
			window.setPixelColour(x, y, color);
		}
	}
}

// Lab 3 Task 5
void drawTextureTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, TextureMap textureMap){
	glm::mat3x3 texture = {{v0.texturePoint.x,v1.texturePoint.x,v2.texturePoint.x},
						   {v0.texturePoint.y,v1.texturePoint.y,v2.texturePoint.y},
							{1,1,1}};
	glm::mat3x3 canvas = {{v0.x,v1.x,v2.x}, {v0.y,v1.y,v2.y},{1,1,1}};
	glm::mat3x3 canvasInverse = glm::inverse(canvas);
	glm::mat3x3 affine = canvasInverse*texture;// should be texture*canvasInverse. why the other way around?
	CanvasTriangle sortedTriangle = sortTriangle(CanvasTriangle(v0,v1,v2));
	CanvasPoint m = findMiddlePoint(sortedTriangle);
	// drawStrokedTriangle(window,sortedTriangle.v0(),m,sortedTriangle.v1(), {255,255,255});
	drawHalfTextureTriangle(window, sortedTriangle.v0(),m,sortedTriangle.v1(),textureMap, affine);
	drawHalfTextureTriangle(window, sortedTriangle.v2(),m,sortedTriangle.v1(),textureMap, affine);
	drawStrokedTriangle(window,sortedTriangle.v0(),sortedTriangle.v1(), sortedTriangle.v2(), {255,255,255});
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

		// Lab 3 task 5
		CanvasPoint v0 = CanvasPoint(160,10);
		CanvasPoint v1 = CanvasPoint(300,230);
		CanvasPoint v2 = CanvasPoint(10,150);
		v0.texturePoint = TexturePoint(195, 5);
		v1.texturePoint = TexturePoint(395, 380);
		v2.texturePoint = TexturePoint(65, 330);
		TextureMap textureMap = TextureMap("texture.ppm");
		drawTextureTriangle(window, v0, v1, v2, textureMap);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
