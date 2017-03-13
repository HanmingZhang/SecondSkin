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

#pragma once

#include <maya/MPxContext.h>
#include <maya/MGlobal.h>
#include <maya/M3dView.h>
#include <maya/MPoint.h>

#include "EDMath.h"

#include <vector>
#include <List>
#include <memory>

class MFnMesh;
class EDAnchor;

class coord {
public:
	short h;
	short v;
	MPoint toMPoint() const;
};

enum EDDrawMode
{
	kDefault,
	kNormal,
	kTangent,
};

struct DrawnCurve
{
	MPoint start;
	MPoint end;
	MString name;
	
	std::shared_ptr<EDAnchor> start_anchor = nullptr;
	std::shared_ptr<EDAnchor> end_anchor = nullptr;

	DrawnCurve(const MPoint & start, const MPoint & end, const MString & name);
};

struct EDAnchor
{
	//TODO: snap to anywhere on current curves
	MPoint point_2D;
	MPoint point_3D;

	EDAnchor() = default;
	EDAnchor(const MPoint & point_2D, const MPoint & point_3D) : point_2D(point_2D), point_3D(point_3D) {}
};

class SecondSkin : public MPxContext
{
public:
	SecondSkin();
	virtual			~SecondSkin();
	void*			creator();

	virtual void toolOnSetup(MEvent & event) override;
	// using Viewport 2.0 of Maya
	virtual MStatus	doPress(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context) override;
	virtual MStatus	doDrag(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context) override;
	virtual MStatus	doRelease(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context) override;
	virtual MStatus drawFeedback(MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context) override;
	virtual void completeAction() override;
	virtual void deleteAction() override;

private:

	void clear_quad_cache();
	void append_stroke(short x, short y);
	void update_anchors();
	void draw_stroke(MHWRender::MUIDrawManager& drawMgr);
	//void draw_anchors(MHWRender::MUIDrawManager& drawMgr);
	bool is_normal(const std::vector<coord> & screen_points, const std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh) const;
	//bool is_tangent() const;
	size_t do_snap(const MPoint & input_end_point);
	MString create_curve(std::vector<coord> & screen_points, MFnMesh* selected_mesh, bool start_known, bool end_known, const MPoint& start_point, MPoint& end_point, bool tangent_mode = false, bool normal_mode = false);
    MString create_surface_from_loop(DrawnCurve& cv);
    void search_loop_from(DrawnCurve &cv, int current_depth, std::list<MString>& loop_list);
    void project_normal(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known);
	void project_tangent(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known);
	void project_contour(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known);
	void project_shell(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known);
	MPoint find_point_nearest_to_mesh(const MFnMesh * selected_mesh, const MPoint & ray_origin, const MVector & ray_direction, const coord & screen_coord, float & ret_height) const;
	void rebuild_kd_2d();
	//void rebuild_kd_3d();
	void rebuild_kd(const MFnMesh * selected_mesh);
    


	bool drawing_quad = true;
	bool firstDraw;
	coord min;
	coord max;
	unsigned maxSize;

	std::vector<coord> stroke;

	//MGlobal::ListAdjustment	listAdjustment;

	M3dView view;
	double normal_threshold = 0.15;
	int tang_samples = 3;
	EDDrawMode drawMode = EDDrawMode::kDefault;

	EDMath::PointCloud<float> mesh_pts;
	EDMath::PointCloud<float> mesh_pts_2d;

	// kd tree for finding nearest point on mesh
	std::unique_ptr<EDMath::KDTree2D> kd_2d = nullptr;

    // TODO: delete these hack
	std::list<MString> prev_curves;
	std::list<std::pair<MPoint, MPoint>> prev_curve_start_end;
    MString prev_surf;

	std::list<DrawnCurve> drawn_curves;

	std::vector<std::shared_ptr<EDAnchor>> anchors; 
	EDMath::PointCloud<float> anchors_2d;
	std::unique_ptr<EDMath::KDTree2D> anchors_kd_2d = nullptr;

	bool first_anchored = false;
    std::shared_ptr<EDAnchor> first_anchor = nullptr;
    std::shared_ptr<EDAnchor> last_anchor = nullptr;
	std::list<MString> drawn_shapes;
};
