#if 0
/**
	This is a simple rendereing engine. 
*/
#endif

#ifndef RENDERER_H
#define RENDERER_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "vec.h"
#include "mat.h" 
#include "utils_math.h"
#include "mesh.h"
#include "scene.h"
#include "camera.h"

typedef enum {
	RP_PERSPECTIVE, 
	RP_ORTHOGRAPHIC
} projection_t;

typedef struct {
	int imgWidth;
	int imgHeight;
	camera_t camera;
	float samplestep;
	cRGB_t * frameBuffer;
	float * zBuffer;
	cRGB_t * texture; //TODO here we need a list of textures...currently we use one as test
	cRGB_t bgcolor;
	projection_t projection;
	float min_z;
	float max_z;
} renderer_t;

typedef struct {
	vec3_t ** pNDC;
	vec3_t ** pRaster;
	float maxx; 
	float maxy; 
	float minx; 
	float miny;
	barycentric_t bc;
	vec3_t pixelSample;
	bool isVisible;
	float over_z;
	float z;
	shape_t * shape;
	renderer_t * renderer;
	int curW;
	int curH;
	unsigned int bufferIndex;
	float cntsamplesuccess;
} renderer_transformation_t;

#if 0
	/**
		render_transformation resource functions.
	*/
#endif
renderer_transformation_t * render_transformation_new();
void render_transformation_free(renderer_transformation_t * rt);

#if 0
	/**
		This functions transforms world to raster. For specific Informations you can see renderer_transformation_t.
	*/
#endif
void process_world_to_raster(renderer_t * renderer, shape_t * shape, renderer_transformation_t * rt);
void process_bounding_box(renderer_transformation_t * rt);
void process_visibility(renderer_transformation_t * rt);
void process_z_calc(renderer_transformation_t * rt);
void process_single_pixel(renderer_transformation_t * rt);
#if 0
	/**
		-render function 
			- worldToCam projection
			- world => NDC conversion
			- NDC => raster conversion
	 
		should add missing vertice or object or something same to parameterlist. Not know yet. prototype data will be set inside.
	 */
#endif
void render_scene(renderer_t * renderer, scene_t * scene);
void render_mesh(renderer_t * renderer, mesh_t * mesh);
void render_shape(renderer_t * renderer, shape_t * shape);

void renderer_clear_frame(renderer_t * renderer);

#if 0
	/**
		Create renderer
	*/
#endif
renderer_t * renderer_new(int imgWidth, int imgHeight, cRGB_t * bgColor);

#if 0
	/**
		Deletes renderer
	*/
#endif
void renderer_free(renderer_t * renderer);

#if 0
	/**
		saves renderer output
	*/
#endif
void renderer_output_ppm(renderer_t * renderer, const char * filename);
#if 0
	/**
		saves z buffer output
	*/
#endif
void renderer_output_z_buffer_ppm(renderer_t * renderer, const char * filename);

#endif