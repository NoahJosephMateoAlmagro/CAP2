#pragma once
/*****************************************************************************/
/* Based on the code written in 2016 by Peter Shirley <ptrshrl@gmail.com>    */
/* Check COPYING.txt for copyright license                                   */
/*****************************************************************************/

#include "Shape.h"

class Sphere : public Shape {
public:
	Sphere(): center(), radius() {}
	Sphere(Vec3 center, float radius) : center(center), radius(radius) {}

	bool collide(const Ray& ray, float t_min, float t_max, CollisionData& cd) const;
	
	Vec3 getCenter() const { return center; }
	float getRadius() const { return radius; }


private:
	Vec3 center;
	float radius;
};
