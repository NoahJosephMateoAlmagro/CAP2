#pragma once

#include "utils.h"
#include "random.h"

#include "Material.h"

class Crystalline : public Material {
public:
	Crystalline(float ri) : ref_idx(ri) {}

	bool scatter(const Ray& r_in, const CollisionData& cd, Vec3& attenuation, Ray& scattered) const;

	float getRefIdx() const { return ref_idx; }

private:
	float ref_idx;
};

