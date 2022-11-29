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
glm::vec3 LIGHT_POINT = {-0.3, 1.6, 1.4}; // Hardcoded

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

std::vector<ModelTriangle> loadOBJFile(std::string OBJfilename, float scalingFactor){
	if (scalingFactor==0.0) throw std::invalid_argument( "Scaling factor shouldn't be zero!" );
	std::vector<ModelTriangle> result;
	std::ifstream OBJFile(OBJfilename);
	std::string line;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> texture_vertices;
	std::vector<Colour> colors;
	std::map<std::string, Colour> mtl;
	std::string MTLfilename = " ";
	int n = 0;
	while (OBJFile.eof()==0){
		n += 1;
		getline(OBJFile, line);
		std::vector<std::string> splitedLine = split(line,' ');
		if(line.empty()) continue;
		else if (splitedLine[0]=="mtllib"){
			// Read .mtl file. Store all values in mtl.
			MTLfilename = splitedLine[1];
			mtl = loadMTLFile(MTLfilename);
		}
		else if (splitedLine[0] == "v"){
			// Load vertices.
			float x = std::stof(splitedLine[1])*scalingFactor;
			float y = std::stof(splitedLine[2])*scalingFactor;
			float z = std::stof(splitedLine[3])*scalingFactor;
			vertices.push_back({x,y,z});
		}
		else if (splitedLine[0] == "vt"){
			// Load Texture points.
			float x = std::stof(splitedLine[1]);
			float y = std::stof(splitedLine[2]);
			texture_vertices.push_back({x,y});

		}
		else if (splitedLine[0] == "f"){	
			// Get 3 points to form a face; get corresponding texture point.
			std::array<glm::vec3, 3> triVer;
			std::array<TexturePoint, 3> text_ver;
			for (int i = 1; i < splitedLine.size(); i++){
				std::vector<std::string> subline = split(splitedLine[i],'/');
				if (subline[1]==""){ // If there's no corresponding texture point.
					float v = std::stof(subline[0]);
					glm::vec3 v0 = vertices[v-1];
					triVer[i-1] = v0;
				} else{ // Get corresponding texture points and store them in text_ver as TexturePoint.
					float v = std::stof(subline[0]);
					float t = std::stof(subline[1]);
					glm::vec3 v0 = vertices[v-1];
					glm::vec2 t0 = texture_vertices[t-1];
					triVer[i-1] = v0;
					text_ver[i-1] = {t0.x,t0.y};
				}
			}
			glm::vec3 v0 = triVer[0];
			glm::vec3 v1 = triVer[1];
			glm::vec3 v2 = triVer[2];
			glm::vec3 v0v1 = v1-v0;
			glm::vec3 v0v2 = v2-v0;
			glm::vec3 normal = glm::normalize(glm::cross(v0v1,v0v2));
			Colour c;
			if (colors.size()!=0) c = colors[colors.size()-1];
			else c = {255,255,255};
			ModelTriangle m;
			m.vertices = triVer;
			m.colour = c;
			m.normal = normal;
			m.texturePoints = text_ver;
			result.push_back(m);
		} 
		// Lab 4 Task 3
		else if (splitedLine[0] == "usemtl" && MTLfilename!=" "){
			std::vector<std::string> splitedLine = split(line,' ');
			Colour c = mtl[splitedLine[1]];
			c.name = splitedLine[1];
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
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj",LOAD_SCALE);
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
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj",LOAD_SCALE);
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
	std::vector<ModelTriangle> modelTriangles = loadOBJFile("cornell-box.obj",LOAD_SCALE);
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
			RayTriangleIntersection rti = RayTriangleIntersection(r, t, modelTriangles[i], i);
			result = rti;
		}
	}
	return result;
}

float getProximity(RayTriangleIntersection rayFromLightPoint){
	glm::vec3 lightToIntersection = rayFromLightPoint.intersectionPoint - LIGHT_POINT;
	float d = glm::length(lightToIntersection);
	float proximity = 20/(4*M_PI*pow(d,2));
	// std::cout << "proximity " <<proximity<< std::endl;
	return glm::clamp(proximity, 0.f, 1.f);
}

float getAngleOfIncidence(glm::vec3 rdFromSurface, glm::vec3 normal){
	float angleOfIncidence = glm::dot(glm::normalize(rdFromSurface),normal);
	// std::cout << "angleOfIncidence " <<glm::clamp(angleOfIncidence, 0.f, 1.f)<< std::endl;
	return glm::clamp(angleOfIncidence, 0.f, 1.f);
}

float getSpecular(glm::vec3 rdFromLightPoint, glm::vec3 rdFromCamera, glm::vec3 normal){
	rdFromCamera = glm::normalize(rdFromCamera);
	glm::vec3 reflectionVector = glm::normalize(rdFromLightPoint - 2*normal*glm::dot(rdFromLightPoint,normal));
	float s = glm::dot(-rdFromCamera,reflectionVector);
	// if (s<0) s = 0.1;
	float specular = pow(s, 65);
	// std::cout << "specular " <<glm::clamp(specular, 0.f, 1.f)<< std::endl;
	return glm::clamp(specular, 0.f, 1.f);
}

std::vector<std::array<glm::vec3,3>> getVertexNormal(std::string OBJFileName, std::vector<ModelTriangle> modelTriangles) {
	std::vector<std::array<glm::vec3,3>> result;
	std::array<std::vector<glm::vec3>,1000> faceNormalsOfVertices;
	std::vector<glm::vec3> faces;
	std::ifstream OBJFile(OBJFileName);
	std::string line;
	int n = 0;
	while(OBJFile.eof()==0){
		getline(OBJFile, line);
		if(line[0]!='f') continue;
		else {
			std::vector<std::string> splitedLine = split(line,' ');
			std::vector<float> vs;
			for (int j = 1; j < splitedLine.size(); j++){
				std::vector<std::string> subline = split(splitedLine[j],'/');
				float v = std::stof(subline[0]);
				vs.push_back(v-1.0);
				glm::vec3 normal = modelTriangles[n].normal;
				float index = v-1.0;
				faceNormalsOfVertices[index].push_back(normal);
			}
			glm::vec3 vss = glm::vec3(vs[0],vs[1],vs[2]);
			faces.push_back(vss);
			n++;
		}
	}
	for (int i = 0; i < modelTriangles.size(); i++){
		std::array<glm::vec3,3> normals;
		for (int j = 0; j<3; j++){
			float f = faces[i][j];
			std::vector<glm::vec3> normalsOfEachV = faceNormalsOfVertices[static_cast< int >(f)];
			glm::vec3 sumOfNormals;
			for (glm::vec3 n : normalsOfEachV) sumOfNormals+=n;
			normals[j] = glm::normalize(sumOfNormals);
		}
		result.push_back(normals);
	}
	return result;
}

// From: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
void Barycentric(glm::vec3 p, ModelTriangle intersectedTriangle, float &u, float &v, float &w) {
	glm::vec3 a = intersectedTriangle.vertices[0];
	glm::vec3 b = intersectedTriangle.vertices[1];
	glm::vec3 c = intersectedTriangle.vertices[2];
    glm::vec3 v0 = b - a;
	glm::vec3 v1 = c - a;
	glm::vec3 v2 = p - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}

glm::vec3 interpolateNormal(glm::vec3 p, RayTriangleIntersection rayFromCamera, std::vector<std::array<glm::vec3,3>> verticesNormals){
	float u;
	float v;
	float w;
	Barycentric(p, rayFromCamera.intersectedTriangle, u, v, w);
	int index = rayFromCamera.triangleIndex;
	std::array<glm::vec3,3> normals = verticesNormals[index];
	glm::vec3 pNormal = u*normals[0]+v*normals[1]+w*normals[2];
	return pNormal;
}

float interpolateBrightness(glm::vec3 p, RayTriangleIntersection rayFromCamera, glm::vec3 brightnesses){
	float u;
	float v;
	float w;
	Barycentric(p, rayFromCamera.intersectedTriangle, u, v, w);
	float brightness = u*brightnesses[0]+v*brightnesses[1]+w*brightnesses[2];
	return glm::clamp(brightness,0.f,1.f);
}

std::vector<glm::vec3> getMultipleLight(glm::vec3 lightPoint, int n_layers){
	std::vector<glm::vec3> result;
	float r = 0.1;
	for (int n = 1; n<n_layers+1; n++){
		int n_points = n*10;
		float theta = 6.28/n_points;
		for (int i = 0; i < n_points; i++){
			float x = lightPoint.x + r*cos(theta);
			float y = lightPoint.y;
			float z = lightPoint.z + r*sin(theta);
			theta+=6.28/n_points;
			result.push_back({x,y,z});
		}
		r+=0.05;
	} 
	return result;
}

void getMirrorColor(RayTriangleIntersection &rti, Colour &c, glm::vec3 intersecP, glm::vec3 reflectionDirection,std::vector<ModelTriangle> modelTriangles){
	while (rti.intersectedTriangle.colour.name == c.name){
		intersecP += rti.distanceFromCamera*reflectionDirection;
		if (rti.distanceFromCamera<1.0e-7) intersecP += 1.0e-6*reflectionDirection;
		rti = getClosestIntersection(intersecP, reflectionDirection, modelTriangles);
	}
	c = rti.intersectedTriangle.colour;
}

float getSoftShadow(RayTriangleIntersection rayFrom, std::vector<ModelTriangle> modelTriangles, float indexC, float brightness_min){
	std::vector<glm::vec3> multipleLights = getMultipleLight(LIGHT_POINT, 5);
	float bri=0.f;
	for (glm::vec3 l : multipleLights){
		glm::vec3 intersectionPoint = rayFrom.intersectionPoint;
		glm::vec3 rdFromLightPoint = intersectionPoint - l;
		RayTriangleIntersection rayFromLightPoint = getClosestIntersection(l, rdFromLightPoint, modelTriangles);
		size_t i = rayFromLightPoint.triangleIndex;
		if (indexC == i) bri +=0.005;
	}
	float brightnessS = glm::clamp(bri, 0.f,brightness_min);
	return brightnessS;
}

void rayTracingRasterisedScene(DrawingWindow &window, std::string OBJFileName){
	std::vector<ModelTriangle> modelTriangles = loadOBJFile(OBJFileName,LOAD_SCALE);
	std::vector<std::array<glm::vec3,3>> verticesNormals = getVertexNormal(OBJFileName,modelTriangles);
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
				float brightness_min = 0.38;
				glm::vec3 normal = interpolateNormal(intersectionPointC,rayFromCamera,verticesNormals);
				// Draw mirror.
				Colour triangleColor = modelTriangles[indexC].colour;
				glm::vec3 reflectionDirection = glm::normalize(rdFromCamera - 2*normal*glm::dot(rdFromCamera, normal));
				RayTriangleIntersection rayFromMirror = getClosestIntersection(intersectionPointC, reflectionDirection, modelTriangles);
				glm::vec3 intersecP = intersectionPointC;
				if (triangleColor.name == "Magenta") getMirrorColor(rayFromMirror,triangleColor,intersecP,reflectionDirection,modelTriangles);
				float brightness = 0.f;
				// Get mirror lighting.
				if (rayFromMirror.intersectionPoint == LIGHT_POINT) brightness +=0.1;
				float proximity = getProximity(rayFromLightPoint);
				float angleOfIncidence = getAngleOfIncidence(-rdFromLightPoint, normal);
				float specular = getSpecular(rdFromLightPoint, rdFromCamera, normal);
				brightness = glm::clamp(proximity*angleOfIncidence+specular, brightness_min,1.f);

				// glm::vec3 intersectionPointM = rayFromMirror.intersectionPoint;
				// glm::vec3 rdFromLightToReflectionPoint = intersectionPointM - LIGHT_POINT;
				// RayTriangleIntersection rayFromLightToReflectionPoint = 
				// 						getClosestIntersection(intersectionPointM, rdFromLightToReflectionPoint, modelTriangles);
				// float indexM = rayFromLightToReflectionPoint.triangleIndex;
				// if (rayFromMirror.triangleIndex != indexM) brightness = getSoftShadow(rayFromMirror,modelTriangles,indexM,brightness_min);
				int r = round(static_cast<float>(triangleColor.red) * brightness);
				int g = round(static_cast<float>(triangleColor.green) * brightness);
				int b = round(static_cast<float>(triangleColor.blue) * brightness);
				uint32_t color = (255 << 24) + (r << 16) + (g << 8) + b;
				if (indexL != indexC){
					float brightnessS = getSoftShadow(rayFromCamera, modelTriangles, indexC, brightness_min);
					Colour c = rayFromCamera.intersectedTriangle.colour;
					if (c.name == "Magenta") getMirrorColor(rayFromMirror,c,intersecP,reflectionDirection,modelTriangles);
					r = round(static_cast<float>(c.red) * brightnessS);
					g = round(static_cast<float>(c.green) * brightnessS);
					b = round(static_cast<float>(c.blue) * brightnessS);
					color = (255 << 24) + (r << 16) + (g << 8) + b;
				}
				window.setPixelColour(j, i, color);
				// Draw shadow.
			}
		}
	}
}

void drawSpherePhong(DrawingWindow &window, std::string OBJFileName){
	std::vector<ModelTriangle> modelTriangles = loadOBJFile(OBJFileName, LOAD_SCALE);
	std::vector<std::array<glm::vec3,3>> verticesNormals = getVertexNormal(OBJFileName, modelTriangles);
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
				float brightness_min = 0.1;
				glm::vec3 normal = interpolateNormal(intersectionPointC, rayFromCamera, verticesNormals);
				// glm::vec3 normal = rayFromLightPoint.intersectedTriangle.normal;
				float proximity = getProximity(rayFromLightPoint);
				float angleOfIncidence = getAngleOfIncidence(-1*rdFromLightPoint, normal);
				float specular = getSpecular(rdFromLightPoint, rdFromCamera, normal);
				// float brightness = proximity;
				// float brightness = angleOfIncidence;
				// float brightness = glm::clamp(specular, brightness_min, 1.f);
				// float brightness = glm::clamp(proximity*angleOfIncidence, brightness_min, 1.f);
				// float brightness = proximity*angleOfIncidence;
				// float brightness = glm::clamp(angleOfIncidence+specular+brightness_min, 0.f, 1.f);
				float brightness = glm::clamp(proximity*angleOfIncidence+specular, brightness_min, 1.f);
				// float brightness = proximity*angleOfIncidence+specular;
				Colour c = {255,0,0};
				int r = round(static_cast<float>(c.red) * brightness);
				uint32_t color = (255 << 24) + (r << 16) + (c.green << 8) + c.blue;
				window.setPixelColour(j, i, color);
			}
		}
	}
}

void drawSphereGouraud(DrawingWindow &window, std::string OBJFileName){
	std::vector<ModelTriangle> modelTriangles = loadOBJFile(OBJFileName, LOAD_SCALE);
	std::vector<std::array<glm::vec3,3>> verticesNormals = getVertexNormal(OBJFileName, modelTriangles);
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
				float brightness_min = 0.1;
				// glm::vec3 normal = interpolateNormal(intersectionPointC, rayFromCamera, verticesNormals);
				
					std::array<glm::vec3,3> normals = verticesNormals[indexL];
					std::vector<float> bs;
					glm::vec3 brightnesses;
					for (int i = 0; i < 3; i++){
						float proximity = getProximity(rayFromLightPoint);
						float angleOfIncidence = getAngleOfIncidence(-1*rdFromLightPoint, normals[i]);
						float specular = getSpecular(rdFromLightPoint, rdFromCamera, normals[i]);
						
						float b = glm::clamp(proximity*angleOfIncidence+specular,brightness_min, 1.f);
						// float b = glm::clamp(proximity+brightness_min, 0.f, 1.f);
						// float b = glm::clamp(proximity*angleOfIncidence, brightness_min, 1.f);
						// float b = glm::clamp(specular, brightness_min, 1.f);
						// float b = glm::clamp(angleOfIncidence, brightness_min, 1.f);
						bs.push_back(b);
					} 
					brightnesses.x = bs[0];
					brightnesses.y = bs[1];
					brightnesses.z = bs[2];
					std::cout<<glm::to_string(brightnesses) << std::endl;
					float brightness = interpolateBrightness(rayFromLightPoint.intersectionPoint, rayFromCamera, brightnesses);
					Colour c = {255,0,0};
					// if (indexC !=indexL) brightness = brightness_min;
					int r = round(static_cast<float>(c.red) * brightness);
					uint32_t color = (255 << 24) + (r << 16) + (c.green << 8) + c.blue;
					window.setPixelColour(j, i, color);
				
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
			rayTracingRasterisedScene(window, "cornell-box.obj");
		}

	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	// loadOBJFile("logo.obj", LOAD_SCALE);
	// loadOBJFile("cornell-box.obj", LOAD_SCALE);
	// loadOBJFile("sphere.obj", LOAD_SCALE);
	cleanBuffer();
	// drawSpherePhong(window, "sphere.obj");
	drawSphereGouraud(window, "sphere.obj");
	// rayTracingRasterisedScene(window, "cornell-box.obj");
	// rayTracingRasterisedScene(window, "logo.obj");
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// rasterisedRender(window);
		// rayTracingRasterisedScene(window, "cornell-box.obj");
		// rayTracingRasterisedScene(window, "logo.obj");
		// drawSpherePhong(window, "sphere.obj");
	
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
