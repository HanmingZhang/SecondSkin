// =============================================================================
//
// EasyDress: a 3D sketching plugin for Maya
// Copyright (C) 2016 
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// =============================================================================

// Created: Mar 29, 2016

#pragma once

#include <nanoflann.hpp>

class MPoint;
class MVector;
class MFnMesh;


namespace EDMath
{
	 MPoint projectOnPlane(const MPoint & point, const MVector & plane_normal
		, const MPoint & ray_origin, const MVector & unit_direction);

	 MVector minimumSkewViewplane(const MVector & ray_direction, const MVector & d);
	 MVector minimumSkewViewplane(const MPoint & camera, const MPoint & p, const MVector & d);
	 double distance_to_mesh(const MFnMesh * selected_mesh, const MPoint & p);

	 template <typename T>
	 struct PointCloud
	 {
		 struct Point
		 {
			 T  x, y, z;
			 Point() = default;
			 Point(T x, T y, T z) : x(x), y(y), z(z) {}
		 };

		 std::vector<Point> pts;

		 inline size_t kdtree_get_point_count() const { return pts.size(); }

		 inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t /*size*/) const
		 {
			 const T d0 = p1[0] - pts[idx_p2].x;
			 const T d1 = p1[1] - pts[idx_p2].y;
			 const T d2 = p1[2] - pts[idx_p2].z;
			 return d0*d0 + d1*d1 + d2*d2;
		 }

		 inline T kdtree_get_pt(const size_t idx, int dim) const
		 {
			 if (dim == 0) return pts[idx].x;
			 else if (dim == 1) return pts[idx].y;
			 else return pts[idx].z;
		 }

		 template <class BBOX>
		 bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

		 void clear()
		 {
			 pts.clear();
		 }
	 };

	 // construct a kd-tree index:
	 typedef nanoflann::KDTreeSingleIndexAdaptor<
		 nanoflann::L2_Simple_Adaptor<float, PointCloud<float> >,
		 PointCloud<float>,
		 3 /* dim */
	 > KDTree3D;

	 // construct a kd-tree index:
	 typedef nanoflann::KDTreeSingleIndexAdaptor<
		 nanoflann::L2_Simple_Adaptor<float, PointCloud<float> >,
		 PointCloud<float>,
		 2 /* dim */
	 > KDTree2D;
}
