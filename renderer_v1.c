#include <stdlib.h>
#include <stdio.h>
#include "renderer_v1.h"

renderer_transformation_t * 
render_transformation_new() {
	return malloc(sizeof(renderer_transformation_t));
}

void 
render_transformation_free(renderer_transformation_t * rt) {
	for(unsigned int ipdnc = 0; ipdnc < rt->shape->cntVertex; ++ipdnc) {
		free(rt->pNDC[ipdnc]);
		free(rt->pRaster[ipdnc]);
	}
	free(rt->pNDC);
	free(rt->pRaster);
	free(rt);
}

void 
process_world_to_raster(renderer_t * renderer, shape_t * shape, renderer_transformation_t * rt) {
	renderer_t * _renderer = renderer;
	shape_t * _shape = shape;
	renderer_transformation_t * _rt = rt;
	camera_t * cam = &_renderer->camera;
	
	_rt->pNDC = malloc(_shape->cntVertex * sizeof(vec3_t));
	_rt->pRaster = malloc(_shape->cntVertex * sizeof(vec3_t));
	
	for(unsigned int ipdnc = 0; ipdnc < _shape->cntVertex; ++ipdnc) {
		_rt->pNDC[ipdnc] = transform_point_new(&cam->transformation, &_shape->vertices[ipdnc]->vec);
		
		_rt->pRaster[ipdnc] = malloc(sizeof(vec3_t));
		
		_rt->pRaster[ipdnc]->x = (_rt->pNDC[ipdnc]->x + 1) * 0.5f * _renderer->imgWidth;
		_rt->pRaster[ipdnc]->y = (1-_rt->pNDC[ipdnc]->y) * 0.5f * _renderer->imgHeight;
		_rt->pRaster[ipdnc]->z = -_rt->pNDC[ipdnc]->z;
	}
	
	if (istriangle(_shape)) {
		_rt->bc.area = place_of_vec3(_rt->pRaster[0], _rt->pRaster[1], _rt->pRaster[2]); 
	}
	
	_rt->renderer = renderer;
	_rt->isVisible = false;
	_rt->over_z = FLT_MAX;
	_rt->shape = shape;
}	

void 
process_bounding_box(renderer_transformation_t * rt) {
	renderer_transformation_t * _rt = rt;
	if(istriangle(_rt->shape)) {
		_rt->maxx = fmin((float)_rt->renderer->imgWidth, fmax(_rt->pRaster[0]->x, fmax(_rt->pRaster[1]->x, _rt->pRaster[2]->x)));
		_rt->maxy = fmin((float)_rt->renderer->imgHeight, fmax(_rt->pRaster[0]->y, fmax(_rt->pRaster[1]->y, _rt->pRaster[2]->y)));
		_rt->minx = fmax(0, fmin(_rt->pRaster[0]->x, fmin(_rt->pRaster[1]->x, _rt->pRaster[2]->x)));
		_rt->miny = fmax(0, fmin(_rt->pRaster[0]->y, fmin(_rt->pRaster[1]->y, _rt->pRaster[2]->y)));
		
	} else if (isline(_rt->shape)) {
		_rt->maxx = fmin((float)_rt->renderer->imgWidth, fmax(_rt->pRaster[0]->x, _rt->pRaster[1]->x));
		_rt->maxy = fmin((float)_rt->renderer->imgHeight, fmax(_rt->pRaster[0]->y, _rt->pRaster[1]->y));
		_rt->minx = fmax(0, fmin(_rt->pRaster[0]->x, _rt->pRaster[1]->x));
		_rt->miny = fmax(0, fmin(_rt->pRaster[0]->y, _rt->pRaster[1]->y));
		
		if ( _rt->minx == _rt->maxx ) {
			_rt->minx -= 2.f; _rt->maxx+=2.f;
		}
		
		if ( _rt->miny == _rt->maxy ) {
			_rt->miny -= 2.f; _rt->maxy+=2.f;
		}
	} else if (ispoint(_rt->shape))
	{
		_rt->maxx = 1.f;
		_rt->maxy = 1.f;
		_rt->minx = 0.f;
		_rt->miny = 0.f;
	}
	
	_rt->curH = _rt->miny;
	_rt->curW = _rt->minx;
	
}

void 
process_visibility(renderer_transformation_t * rt)
{
	renderer_transformation_t * _rt = rt;
	
	_rt->pixelSample.x = _rt->curW + 0.5f;
	_rt->pixelSample.y = _rt->curH + 0.5f;
	_rt->pixelSample.z = 0.f;
	_rt->cntsamplesuccess = 0.f;
	_rt->isVisible = false;
	if (istriangle(_rt->shape) || isline(_rt->shape))
	{
		if (istriangle(_rt->shape))
		{
			is_inside_triangle( _rt->pRaster[0], _rt->pRaster[1], _rt->pRaster[2], &_rt->pixelSample, &_rt->bc);
		} else if(isline(_rt->shape))
		{
			_rt->bc.bc2 = (_rt->pixelSample.x - _rt->pRaster[0]->x ) / (_rt->pRaster[1]->x - _rt->pRaster[0]->x);
		}
		
		//EO is inside after first run test against top left rule
		//if (_rt->bc.inside)
		//{
		
			float maxy = _rt->curH + 1.f;
			float maxx = _rt->curW + 1.f;
			
			float stepstart = 0.5f / _rt->renderer->samplestep; //for st = 2 step is .25  for st = 4 0.125
			float step = 2.f*stepstart; //distance between 
			
			vec3_t subpixel = {0.f, 0.f, 0.f};
			barycentric_t subbc, bcsum;
			bool first = true;
			//printf("stepstart: %f, step %f, curW: %i, curH: %i, bufferindex: %u maxbufferindex: 262143\n", stepstart, step, rt->curW,rt->curH,rt->bufferIndex);
			vec3_t *limitvec;
			subbc.area = _rt->bc.area;
			bcsum.area = _rt->bc.area;
			for (float cury = _rt->curH + stepstart; cury < maxy && cury < _rt->maxy ; cury += step)
			{
				for (float curx = _rt->curW + stepstart; curx < maxx && curx < _rt->maxx ; curx += step)
				{
					subpixel.x = curx;
					subpixel.y = cury;
					//printf("subpx: %f %f\n", curx, cury);
					if ( isline(_rt->shape) )
					{
						float edge = place_of_vec3(_rt->pRaster[0], _rt->pRaster[1], &subpixel);
						limitvec = vec3_sub_new(_rt->pRaster[0], _rt->pRaster[1]);
						vec3_length(limitvec);
						float limit = vec3_length(limitvec) *0.5f;
						if ( ( (edge <= limit && edge >= 0.f) || (edge >= -limit && edge <= 0.f) ) ) 
						{		 
							subbc.inside = true;
						}
						free(limitvec);
					} else
					{
						is_inside_triangle( _rt->pRaster[0], _rt->pRaster[1], _rt->pRaster[2], &subpixel, &subbc);
					}
					
					if(subbc.inside)
					{
						if (first)
						{
							bcsum.bc0 = subbc.bc0;
							bcsum.bc1 = subbc.bc1;
							bcsum.bc2 = subbc.bc2;
							bcsum.w0_12 = subbc.w0_12;
							bcsum.w1_20 = subbc.w1_20;
							bcsum.w2_01 = subbc.w2_01;
							first = false;
							_rt->isVisible = true;
						} else 
						{
							bcsum.bc0 += subbc.bc0;
							bcsum.bc1 += subbc.bc1;
							bcsum.bc2 += subbc.bc2;
							bcsum.w0_12 += subbc.w0_12;
							bcsum.w1_20 += subbc.w1_20;
							bcsum.w2_01 += subbc.w2_01;
						}
						
						++_rt->cntsamplesuccess;
					}
				}
			}

			if(_rt->isVisible)
			{
				_rt->bc.bc0 = bcsum.bc0 / _rt->cntsamplesuccess;
				_rt->bc.bc1 = bcsum.bc1 / _rt->cntsamplesuccess;
				_rt->bc.bc2 = bcsum.bc2 / _rt->cntsamplesuccess;
				_rt->bc.w0_12 = bcsum.w0_12 / _rt->cntsamplesuccess;
				_rt->bc.w1_20 = bcsum.w1_20 / _rt->cntsamplesuccess;
				_rt->bc.w2_01 = bcsum.w2_01 / _rt->cntsamplesuccess;
			}
			//printf("succ: %i\n", _rt->cntsamplesuccess);
		//}//EO is inside after first run test against top left rule
	} else if (ispoint(_rt->shape))
	{
		_rt->isVisible = (_rt->pRaster[0]->x >= 0) &&
					 (_rt->pRaster[0]->x <= _rt->renderer->imgWidth) &&
					 (_rt->pRaster[0]->y >= 0) &&
					 (_rt->pRaster[0]->y <= _rt->renderer->imgHeight);
		if(_rt->isVisible) {
			_rt->curW = _rt->pRaster[0]->x;
			_rt->curH = _rt->pRaster[0]->y;
		}
	}
}

void 
process_z_calc(renderer_transformation_t * rt)
{
	renderer_transformation_t * _rt = rt;
	if (istriangle(_rt->shape))
	{
		_rt->over_z = _rt->pRaster[0]->z * _rt->bc.bc0 +  _rt->pRaster[1]->z * _rt->bc.bc1 + _rt->pRaster[2]->z * _rt->bc.bc2;
	} else if (isline(_rt->shape))
	{
		_rt->over_z = _rt->pRaster[0]->z * (1.f - _rt->bc.bc2) +  _rt->pRaster[1]->z * _rt->bc.bc2;
	} else if (ispoint(_rt->shape))
	{
		_rt->over_z = _rt->pRaster[0]->z;
	}
	
	_rt->z = 1.f / _rt->over_z;
}
//the ORIGINAL ONE DO NOT DELETE :P
void process_texture_(renderer_transformation_t * rt, cRGB_t * curCol)
{
	cRGB_t * curCol1 = curCol;
	renderer_transformation_t * _rt = rt;
	vec2_t tex1, tex2, tex3;
	vec2_copy_dest(&tex1, &_rt->shape->vertices[0]->texCoord);
	vec2_copy_dest(&tex2, &_rt->shape->vertices[1]->texCoord);
	vec2_copy_dest(&tex3, &_rt->shape->vertices[2]->texCoord);
	
	#ifdef debug 
		printf("text coords raw:\n");
		vec2_print(&tex1);
		vec2_print(&tex2);
		vec2_print(&tex3);
	#endif
	
	vec2_mul(&tex1, 1.f/_rt->pRaster[0]->z);
	vec2_mul(&tex2, 1.f/_rt->pRaster[1]->z);
	vec2_mul(&tex3, 1.f/_rt->pRaster[2]->z);
	
	#ifdef debug 
		printf("text coords after div. raster z:\n");
		printf("z: %f div: %f ", _rt->pRaster[0]->z,1.f/_rt->pRaster[0]->z);
		vec2_print(&tex1);
		printf("z: %f div: %f ", _rt->pRaster[1]->z,1.f/_rt->pRaster[1]->z);
		vec2_print(&tex2);
		printf("z: %f div: %f ", _rt->pRaster[2]->z,1.f/_rt->pRaster[2]->z);
		vec2_print(&tex3);
	#endif
	
	float texx = _rt->bc.bc0*tex1.x + 
				 _rt->bc.bc1*tex2.x +
				 _rt->bc.bc2*tex3.x;
	float texy = _rt->bc.bc0*tex1.y + 
				 _rt->bc.bc1*tex2.y +
				 _rt->bc.bc2*tex3.y;
	
	#ifdef debug
		printf("tex coords after mult and add with bc:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	//this cuts texture in perspective...may be right should tested
	texx = (texx / _rt->z)/_rt->z;
	texy = (texy / _rt->z)/_rt->z;
	//can cause number over 1.08f seems a magic number otherwise we must iterate twice for interpolation
	

	#ifdef debug
		printf("tex coords after (texx / _rt->z)/_rt->z:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	texx = interpolate_lin(texx, 0.f, 0.f, 1.08f, 512.f);
	texy = interpolate_lin(texy, 0.f, 0.f, 1.08f, 512.f);
	
	#ifdef debug
		printf("tex coords after miear interpolation:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	int texIndex = floor(texy) * _rt->renderer->imgWidth + floor(texx);
	crgb_crgb_copy(curCol1, &_rt->renderer->texture[texIndex]);
}

void process_texture(renderer_transformation_t * rt, cRGB_t * curCol)
{
	cRGB_t * curCol1 = curCol;
	renderer_transformation_t * _rt = rt;
	
	vec2_t tex1, tex2, tex3;
	vec2_copy_dest(&tex1, &_rt->shape->vertices[0]->texCoord);
	vec2_copy_dest(&tex2, &_rt->shape->vertices[1]->texCoord);
	vec2_copy_dest(&tex3, &_rt->shape->vertices[2]->texCoord);
	
	#ifdef debug 
		printf("text coords raw:\n");
		vec2_print(&tex1);
		vec2_print(&tex2);
		vec2_print(&tex3);
	#endif
	
	vec2_mul(&tex1, 1.f/_rt->pRaster[0]->z);
	vec2_mul(&tex2, 1.f/_rt->pRaster[1]->z);
	vec2_mul(&tex3, 1.f/_rt->pRaster[2]->z);
	
	#ifdef debug 
		printf("text coords after div. raster z:\n");
		printf("z: %f div: %f ", _rt->pRaster[0]->z,1.f/_rt->pRaster[0]->z);
		vec2_print(&tex1);
		printf("z: %f div: %f ", _rt->pRaster[1]->z,1.f/_rt->pRaster[1]->z);
		vec2_print(&tex2);
		printf("z: %f div: %f ", _rt->pRaster[2]->z,1.f/_rt->pRaster[2]->z);
		vec2_print(&tex3);
	#endif
	
	float texx = _rt->bc.bc0*tex1.x + 
				 _rt->bc.bc1*tex2.x +
				 _rt->bc.bc2*tex3.x;
	float texy = _rt->bc.bc0*tex1.y + 
				 _rt->bc.bc1*tex2.y +
				 _rt->bc.bc2*tex3.y;
	
	#ifdef debug
		printf("tex coords after mult and add with bc:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	//this cuts texture in perspective...may be right should tested
	texx = (texx / _rt->z)/_rt->z;
	texy = (texy / _rt->z)/_rt->z;
	//can cause number over 1.08f seems a magic number otherwise we must iterate twice for interpolation
	
    
	#ifdef debug
		printf("tex coords after (texx / _rt->z)/_rt->z:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	texx = interpolate_lin(texx, 0.f, 0.f, 1.08f, 511.f);
	texy = interpolate_lin(texy, 0.f, 0.f, 1.08f, 511.f);
	
	#ifdef debug
		printf("tex coords after miear interpolation:\n");
		printf("texx, texxy:\t%f %f\n", texx, texy);
	#endif
	
	int texIndex = floor(texy) * _rt->renderer->imgWidth + floor(texx);
	crgb_crgb_copy(curCol1, &_rt->renderer->texture[texIndex]);
}

void 
process_color_or_texturing(renderer_transformation_t * rt, unsigned int bi)
{
	renderer_transformation_t * _rt = rt;
	cRGB_t curCol1;
	if (istriangle(_rt->shape))
	{
		if (_rt->shape->texId == -1) 
		{
			crgb_crgb_copy(&curCol1, &_rt->shape->vertices[0]->color);
			cRGB_t curCol2 = { _rt->shape->vertices[1]->color.r, _rt->shape->vertices[1]->color.g , _rt->shape->vertices[1]->color.b };
			cRGB_t curCol3 = { _rt->shape->vertices[2]->color.r, _rt->shape->vertices[2]->color.g , _rt->shape->vertices[2]->color.b };
			
			if ( _rt->renderer->projection == RP_PERSPECTIVE )
			{
				crgb_mul(&curCol1, _rt->bc.bc0);
				crgb_mul(&curCol2, _rt->bc.bc1);
				crgb_mul(&curCol3, _rt->bc.bc2);
			} else {
				crgb_mul(&curCol1, (1.f/_rt->pRaster[0]->z * _rt->bc.bc0));
				crgb_mul(&curCol2, (1.f/_rt->pRaster[1]->z * _rt->bc.bc1));
				crgb_mul(&curCol3, (1.f/_rt->pRaster[2]->z * _rt->bc.bc2));
			}
			
			crgb_crgb_add(&curCol1, &curCol2);
			crgb_crgb_add(&curCol1, &curCol3);
		} else
		{							
			process_texture(_rt, &curCol1);
		}
		
	} else if (isline(_rt->shape))
	{
		crgb_crgb_copy(&curCol1, &_rt->shape->vertices[0]->color);
	} else if (ispoint(_rt->shape))
	{
		crgb_crgb_copy(&curCol1, &_rt->shape->vertices[0]->color);
	}
	
	if (_rt->renderer->samplestep > 1.f )
	{
		float st = _rt->renderer->samplestep*_rt->renderer->samplestep;
		float fail = st - _rt->cntsamplesuccess;
		
		cRGB_t subcolsucc = { curCol1.r * _rt->cntsamplesuccess,
							  curCol1.g * _rt->cntsamplesuccess,
							  curCol1.b * _rt->cntsamplesuccess };
		cRGB_t bufcol = _rt->renderer->frameBuffer[bi/*_rt->bufferIndex*/];
		cRGB_t subcolfail = { bufcol.r * fail,
							  bufcol.g * fail,
							  bufcol.b * fail};
		crgb_crgb_add(&subcolsucc, &subcolfail);
		crgb_mul(&subcolsucc, 1.f/st);
		crgb_crgb_copy(&curCol1, &subcolsucc);
	}
	
	crgb_crgb_copy(&_rt->renderer->frameBuffer[bi/*_rt->bufferIndex*/], &curCol1);
}

void process_single_pixel(renderer_transformation_t * rt)
{
	renderer_transformation_t * _rt = rt;
	process_visibility(_rt);
	if(_rt->isVisible) {

		process_z_calc(_rt);
		
		//_rt->bufferIndex = _rt->curH * _rt->renderer->imgWidth + _rt->curW;
		unsigned int bi = _rt->curH * _rt->renderer->imgWidth + _rt->curW;
		if ( _rt->z < _rt->renderer->zBuffer[bi/*_rt->bufferIndex*/] )
		{
			_rt->renderer->zBuffer[bi/*_rt->bufferIndex*/] = _rt->z;
		
			_rt->renderer->min_z = fminf(_rt->renderer->min_z, _rt->z);
			_rt->renderer->max_z = fmaxf(_rt->renderer->max_z, _rt->z);

			process_color_or_texturing(_rt, bi);
		}
	}
}

void render_shape_raw(renderer_t * renderer, shape_t * shape)
{
	renderer_transformation_t * rt = render_transformation_new();
	process_world_to_raster(renderer, shape, rt);
	process_bounding_box(rt);
	
	for(; rt->curH < rt->maxy; ++rt->curH)
	{
		for(rt->curW = rt->minx; rt->curW < rt->maxx; ++rt->curW)
		{
			process_single_pixel(rt);
		}
	}
	
	render_transformation_free(rt);
}



void render_mesh_raw(renderer_t * renderer, mesh_t * mesh)
{
	for(unsigned int cntShape = 0; cntShape < mesh->cntShapes ; ++cntShape) {
		render_shape_raw(renderer, mesh->shapes[cntShape]);
	} 
}

void 
render_scene_raw(renderer_t * renderer, scene_t * scene)
{	
	for(int cntMesh = 0; cntMesh < scene->cntMesh ; ++cntMesh) {
		render_mesh_raw(renderer, scene->meshes[cntMesh]);
	}
}

void render_shape(renderer_t * renderer, shape_t * shape)
{
	float samples = renderer->samplestep;
	renderer->samplestep = 1.f;
	render_shape_raw(renderer, shape);
	if (samples > renderer->samplestep) {
		renderer->samplestep = samples;
		render_shape_raw(renderer, shape);
	}
}

void render_mesh(renderer_t * renderer, mesh_t * mesh)
{
	float samples = renderer->samplestep;
	renderer->samplestep = 1.f;
	render_mesh_raw(renderer, mesh);
	if (samples > renderer->samplestep) {
		renderer->samplestep = samples;
		render_mesh_raw(renderer, mesh);
	}
}

void 
render_scene(renderer_t * renderer, scene_t * scene)
{	
	float samples = renderer->samplestep;
	renderer->samplestep = 1.f;
	render_scene_raw(renderer, scene);
	if (samples > renderer->samplestep) {
		renderer->samplestep = samples;
		render_scene_raw(renderer, scene);
	}
}

void renderer_clear_frame(renderer_t * renderer)
{
	int buffersize = renderer->imgWidth * renderer->imgHeight;
	for(int i = 0; i < buffersize; ++i)
	{
		crgb_crgb_copy(&renderer->frameBuffer[i], &renderer->bgcolor);
	}
	
	for(int i = 0; i < buffersize; ++i)
	{
		renderer->zBuffer[i] = FLT_MAX;
	}
	renderer->min_z = FLT_MAX;
	renderer->max_z = 0.f;
}

renderer_t * 
renderer_new(int imgWidth, int imgHeight, cRGB_t * bgColor)
{
	renderer_t * newrenderer = malloc(sizeof(renderer_t));
	newrenderer->projection = RP_ORTHOGRAPHIC;
	newrenderer->texture = NULL;
	newrenderer->imgWidth = imgWidth;
	newrenderer->imgHeight = imgHeight;
	int buffersize = imgWidth * imgHeight;
	newrenderer->frameBuffer = malloc(buffersize * sizeof(cRGB_t));
	newrenderer->zBuffer = malloc(buffersize * sizeof(float));
	crgb_crgb_copy(&newrenderer->bgcolor, bgColor);
	renderer_clear_frame(newrenderer);
	newrenderer->min_z = FLT_MAX;
	newrenderer->max_z = 0.f;
	newrenderer->samplestep = 1.f;
	return newrenderer;
}

void 
renderer_free(renderer_t * renderer)
{
	free(renderer->frameBuffer);
	free(renderer->zBuffer);
	if (renderer->texture != NULL) {
		free(renderer->texture);
	}
}


void 
renderer_output_ppm(renderer_t * renderer, const char * filename)
{
	int i, j;
    FILE *fp = fopen(filename, "wb"); /* b - binary mode */
    (void) fprintf(fp, "P6\n%d %d\n255\n", renderer->imgWidth, renderer->imgHeight);
	
    for (j = 0; j < renderer->imgHeight; ++j)
    {
	  for (i = 0; i < renderer->imgWidth; ++i)
	  {
		static unsigned char color[3];
	
		color[0] = (unsigned char)interpolate_lin(renderer->frameBuffer[j * renderer->imgWidth + i].r, 0.f, 0.f, 1.1f, 255.f);
		color[1] = (unsigned char)interpolate_lin(renderer->frameBuffer[j * renderer->imgWidth + i].g, 0.f, 0.f, 1.1f, 255.f);
		color[2] = (unsigned char)interpolate_lin(renderer->frameBuffer[j * renderer->imgWidth + i].b, 0.f, 0.f, 1.1f, 255.f);
		(void) fwrite(color, 1, 3, fp);
	  }
    }
    (void) fclose(fp);
}

void 
renderer_output_z_buffer_ppm(renderer_t * renderer, const char * filename)
{
	int i, j;
    FILE *fp = fopen(filename, "wb"); /* b - binary mode */
    (void) fprintf(fp, "P6\n%d %d\n255\n", renderer->imgWidth, renderer->imgHeight);

	unsigned char curcol;
	float curz;
    for (j = 0; j < renderer->imgHeight; ++j)
    {
	  for (i = 0; i < renderer->imgWidth; ++i)
	  {
		static unsigned char color[3];
		
		curz = renderer->zBuffer[j * renderer->imgWidth + i];
		if ( curz != FLT_MAX )
		{
			curcol = (unsigned char)interpolate_lin(curz, renderer->max_z, 0.f, renderer->min_z, 255.f);
		} else {
			curcol = 0.f;
		}
		
		color[0] = curcol;  
		color[1] = curcol;  
		color[2] = curcol;
		(void) fwrite(color, 1, 3, fp);
	  }
    }
    (void) fclose(fp);
}
