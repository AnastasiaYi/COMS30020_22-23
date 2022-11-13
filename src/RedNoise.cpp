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
// Lab 4 Task 2
#include <ModelTriangle.h>
#include <iostream>
// Lab 4 Task 3
#include <map>
// Lab 6 Task 2
#include <RayTriangleIntersection.h>

// #define WIDTH 320
// #define HEIGHT 240
// #define WIDTH 640
// #define HEIGHT 480
#define WIDTH 960
#define HEIGHT 720

// Lab 4
#define LOAD_SCALE 0.35
#define IMAGE_SCALE 300
#define FOCAL_LENGTH 2.0
float DEPTH_BUFFER[HEIGHT][WIDTH];

// Lab 5
glm::vec3 CAMERA_POSITION = glm::vec3(0,0,4);
glm::mat3 CAMERA_ORIENTATION = glm::mat3(1,0,0,0,1,0,0,0,1);
float cosine = cos(0.01);
float sine = sin(0.01);

// Lab 6
glm::vec3 LIGHT_POINT = {0.0, 0.8,0.0}; // Hardcoded

void cleanBuffer() {
	for(int i = 0; i < HEIGHT; i++) {
		for(int j = 0; j < WIDTH; j++) {
			DEPTH_BUFFER[i][j] = 0;
		}
	}
}

std::vector<glm::vec2> interpolateTwoElementValues(glm::vec2 from, glm::vec2 to, float numberOfValues){
	std::vector<glm::vec2> result;
	glm::vec2 stepSize = (to-from)/(numberOfValues);
	for (float i = 0; i < numberOfValues; i++){
		result.push_back(from + i*stepSize);
	}
	return result;
}

CanvasTriangle generateRandomTriangle(DrawingWindow &window){
	CanvasPoint v0 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasPoint v1 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasPoint v2 = {float(rand()%window.width), float(rand()%window.height)};
	CanvasTriangle aLovelyTriangle = CanvasTriangle(v0, v1, v2);
	return aLovelyTriangle;
}

std::vector<float> interpolateSingleFloats(float from, float to, float numberOfValues) {
	std::vector<float> result;
	float stepSize = (to - from)/(numberOfValues-1);
	for (int i = 0; i < numberOfValues; i ++) {
		result.push_back(from + i*stepSize);
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, float numberOfValues){
	std::vector<glm::vec3> result;
	glm::vec3 stepSize = (to - from)/(numberOfValues);
	for (float i = 0; i < numberOfValues; i++){ // i must be float. Can't multiply int with float in c++.
		result.push_back(from + i*stepSize);
	}
	return result;
}

// TODO: fix skipped line.
void drawLine (DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour color){
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = 1/to.depth - 1/from.depth;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff))+5;
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	float zStepSize = zDiff/numberOfSteps;
	for (float i = 0.0; i < numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		float depth = (1/from.depth + (zStepSize*i));
		if ((DEPTH_BUFFER[int(floor(y))][int(floor(x))] <= depth)|(DEPTH_BUFFER[int(y)][int(x)] == 0)){
			uint32_t c = (255 << 24) + (color.red << 16) + (color.green << 8) + color.blue;
			window.setPixelColour(int(floor(x)), int(floor(y)), c);
			DEPTH_BUFFER[int(floor(y))][int(floor(x))] = depth;
		}
	}
}

void drawStrokedTriangle (DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, Colour color){
	drawLine(window, v0, v1, color);
	drawLine(window, v1, v2, color);
	drawLine(window, v0, v2, color);
}

CanvasTriangle sortTriangle (CanvasTriangle triangle){
	while (!((triangle.v0().y<=triangle.v1().y)&&(triangle.v1().y<=triangle.v2().y))){
		if (triangle.v0().y>=triangle.v1().y) std::swap(triangle.v0(),triangle.v1());
		if (triangle.v1().y>=triangle.v2().y) std::swap(triangle.v1(), triangle.v2());
	}
	CanvasTriangle sortedTriangle = {triangle.v0(), triangle.v1(), triangle.v2()};
	return sortedTriangle;
}

CanvasPoint findMiddlePoint (CanvasTriangle triangle){
	float x0 = triangle.v0().x;
	float x2 = triangle.v2().x;
	float z0 = 1/triangle.v0().depth;
	float z2 = 1/triangle.v2().depth;
	float y0 = triangle.v0().y;
	float y1 = triangle.v1().y;
	float y2 = triangle.v2().y;
	float x = x2+(y1-y2)*(x0-x2)/(y0-y2);
	float d = z0-(y0-y1)*(z0-z2)/(y0-y2);
	CanvasPoint c = CanvasPoint(x, y1, 1/d);
	return c;
}

void drawHalfTriangle (DrawingWindow &window, CanvasPoint h0, CanvasPoint h1, CanvasPoint h2, Colour color){
	glm::vec3 from = {h0.x, h0.y, 1/h0.depth};
	glm::vec3 to1 = {h1.x, h1.y, 1/h1.depth};
	glm::vec3 to2 = {h2.x, h1.y,1/h2.depth};
	float numberOfValues = abs(h1.y - h0.y)+5;
	std::vector<glm::vec3> h0h1  = interpolateThreeElementValues(from, to1, numberOfValues);
	std::vector<glm::vec3> h0h2  = interpolateThreeElementValues(from, to2, numberOfValues);
	for (int i = 0; i < numberOfValues; i++){
		CanvasPoint t = CanvasPoint(h0h2[i][0], h0h2[i][1], 1/h0h2[i][2]);
		CanvasPoint f = CanvasPoint(h0h1[i][0], h0h1[i][1], 1/h0h1[i][2]);
		drawLine(window, f,t, color);
	}
}

void drawFilledTriangle (DrawingWindow &window, CanvasTriangle triangle, Colour color){
	// Colour white = {255, 255, 255};
	CanvasTriangle sortedTriangle = sortTriangle(triangle);
	CanvasPoint m = findMiddlePoint(sortedTriangle);
	drawHalfTriangle(window, sortedTriangle.v0(), m, sortedTriangle.v1(), color);
	drawHalfTriangle(window, sortedTriangle.v2(), m, sortedTriangle.v1(), color);
	drawStrokedTriangle(window, sortedTriangle.v0(), sortedTriangle.v1(), sortedTriangle.v2(), color);
}

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

void drawTextureTriangle(DrawingWindow &window, CanvasPoint v0, CanvasPoint v1, CanvasPoint v2, TextureMap textureMap){
	glm::mat3x3 texture = {{v0.texturePoint.x,v1.texturePoint.x,v2.texturePoint.x},
						   {v0.texturePoint.y,v1.texturePoint.y,v2.texturePoint.y},
							{1,1,1}};
	glm::mat3x3 canvas = {{v0.x,v1.x,v2.x}, {v0.y,v1.y,v2.y},{1,1,1}};
	glm::mat3x3 canvasInverse = glm::inverse(canvas);
	glm::mat3x3 affine = canvasInverse*texture;// should be texture*canvasInverse. why the other way around?
	CanvasTriangle sortedTriangle = sortTriangle(CanvasTriangle(v0,v1,v2));
	CanvasPoint m = findMiddlePoint(sortedTriangle);
	drawHalfTextureTriangle(window, sortedTriangle.v0(),m,sortedTriangle.v1(),textureMap, affine);
	drawHalfTextureTriangle(window, sortedTriangle.v2(),m,sortedTriangle.v1(),textureMap, affine);
	drawStrokedTriangle(window,sortedTriangle.v0(),sortedTriangle.v1(), sortedTriangle.v2(), {255,255,255});
}

std::map<std::string, Colour> loadMTLFile(std::string filename) {
	std::map<std::string, Colour> result;
	std::ifstream MTLFile(filename);
	std::string line;
	std::vector<std::string> name;
	while (MTLFile.eof()==0){
		getline(MTLFile, line);
		if (line.empty()) continue;
		std::vector<std::string> splitedLine = split(line,' ');
		if (splitedLine[0]=="newmtl") name.push_back(splitedLine[1]);
		else if (splitedLine[0]=="Kd"){
			float r = std::stof(splitedLine[1])*255;
			float g = std::stof(splitedLine[2])*255;
			float b = std::stof(splitedLine[3])*255;
			Colour c = {name[name.size()-1],static_cast<int>(r),static_cast<int>(g),static_cast<int>(b)};
			result[c.name] = c;
		}
	}
	return result;
}

std::vector<ModelTriangle> loadOBJFile(std::string OBJfilename, std::string MTLfilename, float scalingFactor){
	if (scalingFactor==0.0) throw std::invalid_argument( "Scaling factor shouldn't be zero!" );
	std::vector<ModelTriangle> result;
	std::ifstream OBJFile(OBJfilename);
	std::string line;
	std::vector<glm::vec3> vertices;
	std::vector<Colour> colors;
	std::map<std::string, Colour> mtl = loadMTLFile(MTLfilename);
	while (OBJFile.eof()==0){
		getline(OBJFile, line);
		if(line.empty()) continue;
		else if (line[0] == 'v'){
			std::vector<std::string> splitedLine = split(line,' ');
			float x = std::stof(splitedLine[1])*scalingFactor;
			float y = std::stof(splitedLine[2])*scalingFactor;
			float z = std::stof(splitedLine[3])*scalingFactor;
			vertices.push_back({x,y,z}); // why do we have to set x,y,z but not pass stof() in directly.
		}
		else if (line[0] == 'f'){	
			line.erase(remove(line.begin(), line.end(), '/'), line.end());
			std::vector<std::string> splitedLine = split(line,' ');
			float x = std::stof(splitedLine[1]);
			float y = std::stof(splitedLine[2]);
			float z = std::stof(splitedLine[3]);
			glm::vec3 v0 = vertices[x-1];
			glm::vec3 v1 = vertices[y-1];
			glm::vec3 v2 = vertices[z-1];
			glm::vec3 v0v1 = v0-v1;
			glm::vec3 v1v2 = v1-v2;
			glm::vec3 normal = glm::cross(v0v1,v1v2);
			ModelTriangle m = ModelTriangle(v0, v1, v2, colors[colors.size()-1]);
			m.normal = normal;
			result.push_back(m);
		} 
		// Lab 4 Task 3
		else if (line[0] == 'u'){
			std::vector<std::string> splitedLine = split(line,' ');
			Colour c = mtl[splitedLine[1]];
			colors.push_back(c);
		}
    }
	return result;
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength, float imagePlaneScaling){
	glm::vec3 vertexP = (vertexPosition-cameraPosition)*CAMERA_ORIENTATION;
	// glm::vec3 vertexP = vertexPosition-cameraPosition;

	float u = -1*imagePlaneScaling*focalLength * vertexP.x/vertexP.z + WIDTH/2;
	float v = imagePlaneScaling*focalLength * vertexP.y/vertexP.z + HEIGHT/2;
	float d = 1/vertexP.z;
	return CanvasPoint(u,v,d);
}

void pointcloud(DrawingWindow &window){
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj","cornell-box.mtl",LOAD_SCALE);
	uint32_t color = (255 << 24) + (255 << 16) + (255 << 8) + 255;
	for (int i = 0; i < modelTriangles.size(); i++){
		std::array<glm::vec3,3> vertices = modelTriangles[i].vertices;
		for (int j = 0; j < vertices.size(); j++){
			CanvasPoint c = getCanvasIntersectionPoint(CAMERA_POSITION, vertices[j], FOCAL_LENGTH, IMAGE_SCALE);
			window.setPixelColour(c.x,c.y,color);
		}
	}
}

void wireframeRender(DrawingWindow &window){
	cleanBuffer();
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj","cornell-box.mtl",LOAD_SCALE);
	for (int i = 0; i < modelTriangles.size(); i++){
		std::array<glm::vec3,3> vertices = modelTriangles[i].vertices;
		glm::vec3 v0 = vertices[0];
		glm::vec3 v1 = vertices[1];
		glm::vec3 v2 = vertices[2];
		CanvasPoint c0 = getCanvasIntersectionPoint(CAMERA_POSITION, v0, FOCAL_LENGTH, IMAGE_SCALE);
		CanvasPoint c1 = getCanvasIntersectionPoint(CAMERA_POSITION, v1, FOCAL_LENGTH, IMAGE_SCALE);
		CanvasPoint c2 = getCanvasIntersectionPoint(CAMERA_POSITION, v2, FOCAL_LENGTH, IMAGE_SCALE);
		drawStrokedTriangle(window,c0,c1,c2,{255,255,255});
	}
}

void rasterisedRender(DrawingWindow &window){
	cleanBuffer();
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj","cornell-box.mtl",LOAD_SCALE);
	for (int i = 0; i < modelTriangles.size(); i++){
		std::array<glm::vec3,3> vertices = modelTriangles[i].vertices;
		glm::vec3 v0 = vertices[0];
		glm::vec3 v1 = vertices[1];
		glm::vec3 v2 = vertices[2];
		CanvasPoint c0 = getCanvasIntersectionPoint(CAMERA_POSITION, v0, FOCAL_LENGTH, IMAGE_SCALE);
		CanvasPoint c1 = getCanvasIntersectionPoint(CAMERA_POSITION, v1, FOCAL_LENGTH, IMAGE_SCALE);
		CanvasPoint c2 = getCanvasIntersectionPoint(CAMERA_POSITION, v2, FOCAL_LENGTH, IMAGE_SCALE);
		CanvasTriangle t = CanvasTriangle(c0,c1,c2);
		drawFilledTriangle(window,t,modelTriangles[i].colour);
	}
}

void lookAt(glm::vec3 lookAtPoint){
	glm::vec3 forward = glm::normalize(CAMERA_POSITION-lookAtPoint);
	glm::vec3 right = glm::cross({0,1,0},forward);
	glm::vec3 up = glm::cross(forward,right);
	CAMERA_ORIENTATION = glm::mat3(right,up,forward);
}

void cameraMoving(DrawingWindow &window, float step, int index){
	window.clearPixels();
	CAMERA_POSITION[index] += step;
	lookAt({0,0,0});
	rasterisedRender(window);
}

void cameraRotation(glm::mat3 m, DrawingWindow &window){
	cleanBuffer();
	window.clearPixels();
	CAMERA_POSITION = m*CAMERA_POSITION;
	lookAt({0,0,0});
	rasterisedRender(window);
}

void orbit(DrawingWindow &window,SDL_Event event){
	glm::mat3 m = glm::mat3(cosine,0, -1*sine, 
								0, 1, 0, 
								sine, 0, cosine);
	while (true){
		if (window.pollForInputEvents(event)) if (event.key.keysym.sym == 'q') break;
		cameraRotation(m, window);
		window.renderFrame();
		
	}
}

RayTriangleIntersection getClosestIntersection(glm::vec3 startingPoint, glm::vec3 rayDirection, std::vector<ModelTriangle> modelTriangles){
	RayTriangleIntersection result;
	result.distanceFromCamera = INFINITY;
	for (int i = 0; i < modelTriangles.size(); i++) {
		glm::vec3 p0 = modelTriangles[i].vertices[0];
		glm::vec3 p1 = modelTriangles[i].vertices[1];
		glm::vec3 p2 = modelTriangles[i].vertices[2];
		glm::vec3 e0 = p1 - p0;
		glm::vec3 e1 = p2 - p0;
		glm::vec3 SPVector = startingPoint - p0;
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 TUVMatrix = glm::inverse(DEMatrix) * SPVector;
		float t = TUVMatrix[0];
		float u = TUVMatrix[1];
		float v = TUVMatrix[2];
		if (t < result.distanceFromCamera && t > 0 && u >= 0 && u <= 1 && v >= 0 && v <= 1 && (u + v) <= 1) {
			glm::vec3 r = startingPoint + t*rayDirection;
			RayTriangleIntersection rti = RayTriangleIntersection(r, TUVMatrix[0], modelTriangles[i], i);
			result = rti;
		}
	}
	return result;
}

float getProximity(RayTriangleIntersection rayFromLightPoint){
	glm::vec3 lightToIntersection = rayFromLightPoint.intersectionPoint - LIGHT_POINT;
	float d = glm::length(lightToIntersection);
	float proximity = 22/(4*M_PI*pow(d,2));
	if (proximity > 1) proximity = 1;
	else if (proximity < 0) proximity = 0;
	return proximity;
}

float getAngleOfIncidence(ModelTriangle intersectedTriangle,glm::vec3 rdFromSurface){
	glm::vec3 normal = intersectedTriangle.normal;
	float normalLength = glm::length(normal);
	float rayLength = glm::length(rdFromSurface);
	float angleOfIncidence = glm::dot(normal, rdFromSurface)/(normalLength*rayLength);
	return angleOfIncidence;
}

float getSpecular(ModelTriangle intersectedTriangle, glm::vec3 rdFromLightPoint, glm::vec3 rdFromCamera){
	glm::vec3 normal = intersectedTriangle.normal;
	float normalLength = glm::length(normal);
	float rayFromLightPointLength = glm::length(rdFromLightPoint);
	float rayFromCameraLength = glm::length(rdFromCamera);
	normal /= normalLength;
	rdFromLightPoint /= rayFromLightPointLength;
	rdFromCamera /= rayFromCameraLength;
	glm::vec3 reflectionVector = rdFromLightPoint - 2*normal*glm::dot(rdFromLightPoint,normal);
	float s = glm::dot(-1*rdFromCamera,reflectionVector);
	// if (s<0) s = 0.1;
	float specular = pow(s, 2);
	// if (specular>1) specular = 1;
	// else if (specular<0) specular = 0;
	// std::cout << specular << std::endl;
	return specular;
}

Colour getColour(ModelTriangle intersectedTriangle, RayTriangleIntersection rayFromLightPoint, glm::vec3 rdFromLightPoint, glm::vec3 rdFromCamera, float brightness_min){
	glm::vec3 rdFromSurface = -1*rdFromLightPoint;
	float proximity = getProximity(rayFromLightPoint);
	float angleOfIncidence = getAngleOfIncidence(intersectedTriangle,rdFromSurface);
	float specular = getSpecular(intersectedTriangle, rdFromLightPoint, rdFromCamera);
	// float brightness = proximity;
	// float brightness = angleOfIncidence;
	// float brightness = specular;
	float brightness = proximity*angleOfIncidence;
	// float brightness = proximity*angleOfIncidence*specular;
	if (brightness>1) brightness = 1;
	else if (brightness<brightness_min) brightness = brightness_min;
	Colour triangleColor = intersectedTriangle.colour;
	int r = round(static_cast<float>(triangleColor.red) * brightness);
	int g = round(static_cast<float>(triangleColor.green) * brightness);
	int b = round(static_cast<float>(triangleColor.blue) * brightness);
	return Colour(r,g,b);
}

void rayTracingRasterisedScene(DrawingWindow &window){
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj","cornell-box.mtl",LOAD_SCALE);
	for (float i = 0; i < HEIGHT; i++){
		for (float j = 0; j < WIDTH; j++){
			float z = CAMERA_POSITION.z - FOCAL_LENGTH;
			float x = (j-WIDTH/2)*z/(FOCAL_LENGTH*IMAGE_SCALE);
			float y = -1*(i-HEIGHT/2)*z/(FOCAL_LENGTH*IMAGE_SCALE);
			glm::vec3 canvasTo3D = {x,y,z};
			glm::vec3 rdFromCamera = canvasTo3D - CAMERA_POSITION;
			RayTriangleIntersection rayFromCamera = getClosestIntersection(CAMERA_POSITION, rdFromCamera, modelTriangles);
			size_t indexC = rayFromCamera.triangleIndex;

			if (rayFromCamera.distanceFromCamera!=INFINITY){
				glm::vec3 intersectionPointC = rayFromCamera.intersectionPoint;
				glm::vec3 rdFromLightPoint = intersectionPointC - LIGHT_POINT;
				RayTriangleIntersection rayFromLightPoint = getClosestIntersection(LIGHT_POINT, rdFromLightPoint, modelTriangles);
				size_t indexL = rayFromLightPoint.triangleIndex;
				float brightness_min = 0.2;
				if (indexL == indexC){
					Colour c = getColour(rayFromCamera.intersectedTriangle, rayFromLightPoint, rdFromLightPoint, rdFromCamera, brightness_min);
					uint32_t color = (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
					window.setPixelColour(j, i, color);
				}else{
					Colour c = rayFromCamera.intersectedTriangle.colour;
					int r = round(static_cast<float>(c.red) * brightness_min);
					int g = round(static_cast<float>(c.green) * brightness_min);
					int b = round(static_cast<float>(c.blue) * brightness_min);
					uint32_t color = (255 << 24) + (r << 16) + (g << 8) + b;
					window.setPixelColour(j, i, color);
				}
			}
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	Colour color = {rand()%255, rand()%255, rand()%255};
	if (event.type == SDL_KEYDOWN) {
		// Lab 5 Task 3 changed to manipulate camera position
		if (event.key.keysym.sym == SDLK_LEFT){
			cameraMoving(window, 0.17, 0);
		} else if (event.key.keysym.sym == SDLK_RIGHT){
			cameraMoving(window, -0.17, 0);
		} else if (event.key.keysym.sym == SDLK_UP){
			cameraMoving(window, 0.17, 1);
		} else if (event.key.keysym.sym == SDLK_DOWN){
			cameraMoving(window, -0.17, 1);
		} else if (event.key.keysym.sym == 'j'){
			cameraMoving(window, 0.17, 2);
		} else if (event.key.keysym.sym == 'l'){
			cameraMoving(window, -0.17, 2);
		}
		// Lab 5 Task 3 rotate
		// Rotation about Y
		else if (event.key.keysym.sym == 'a'){
			// rotate to right
			glm::mat3 m = glm::mat3(cosine,0, -1*sine, 
									0, 1, 0, 
									sine, 0, cosine);
			cameraRotation(m,window);
		} else if (event.key.keysym.sym == 'd'){
			// rotate to left
			glm::mat3 m = glm::mat3(cosine,0, sine, 
									0, 1, 0, 
									-1*sine, 0, cosine);
			cameraRotation(m,window);		
		}
		// Rotation about X
		else if (event.key.keysym.sym == 'w'){
			// rotate downwards
			glm::mat3 m = glm::mat3(1,0, 0, 
									0, cosine, sine, 
									0, -1*sine, cosine);
			cameraRotation(m,window);
		} else if (event.key.keysym.sym == 's'){
			// rotate upwards
			glm::mat3 m = glm::mat3(1,0, 0, 
									0, cosine, -1*sine, 
									0, sine, cosine);
			cameraRotation(m,window);
		} else if (event.key.keysym.sym == 'r'){
			orbit(window,event);
		}
		// Lab 3 Task 3
		else if (event.key.keysym.sym == 'u') {
			Colour color = {rand()%255, rand()%255, rand()%255};
			CanvasTriangle triangle = generateRandomTriangle(window);
			drawStrokedTriangle(window, triangle.v0(),triangle.v1(),triangle.v2(), color);
		} else if (event.key.keysym.sym == 'f') {
			CanvasTriangle triangle = generateRandomTriangle(window);
			drawFilledTriangle(window, triangle, color);
		} else if (event.key.keysym.sym == 't') {
			// Lab 3 Task 5/6
			CanvasPoint v0 = CanvasPoint(160,10);
			CanvasPoint v1 = CanvasPoint(300,230);
			CanvasPoint v2 = CanvasPoint(10,150);
			v0.texturePoint = TexturePoint(195, 5);
			v1.texturePoint = TexturePoint(395, 380);
			v2.texturePoint = TexturePoint(65, 330);
			TextureMap textureMap = TextureMap("texture.ppm");
			drawTextureTriangle(window, v0, v1, v2, textureMap);
		}
		// Lab 6 Task 7
		else if (event.key.keysym.sym == 'z'){
			window.clearPixels();
			wireframeRender(window);
		}else if (event.key.keysym.sym == 'x'){
			window.clearPixels();
			rasterisedRender(window);
		}else if (event.key.keysym.sym == 'c'){
			window.clearPixels();
			rayTracingRasterisedScene(window);
		}

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	cleanBuffer();
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// rasterisedRender(window);
		rayTracingRasterisedScene(window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
