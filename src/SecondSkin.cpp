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

// The main tool for sketching in viewport.
//
// Created: Mar 4, 2016


#include "SecondSkin.h"

#include <maya/MItSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MDagPath.h>
#include <maya/MFnTransform.h>
#include <maya/MUIDrawManager.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnMesh.h>
#include <maya/MPointArray.h>

#include <nanoflann.hpp>
#include "EDMath.h"

#include <string>
#include <list>
#include <vector>

DrawnCurve::DrawnCurve(const MPoint & start, const MPoint & end, const MString & name)
	: start(start), end(end), name(name)
{}


const int initialSize = 1024;
const int increment = 256;
const char helpString[] = "drag mouse to draw strokes";
extern "C" int xycompare(coord *p1, coord *p2);
int xycompare(coord *p1, coord *p2)
{
	if (p1->v > p2->v)
		return 1;
	if (p2->v > p1->v)
		return -1;
	if (p1->h > p2->h)
		return 1;
	if (p2->h > p1->h)
		return -1;

	return 0;
}

SecondSkin::SecondSkin()
{
	setTitleString("SecondSkin");
}

SecondSkin::~SecondSkin() {}

void* SecondSkin::creator()
{
	return new SecondSkin;
}

void SecondSkin::clear_quad_cache()
{
	prev_curves.clear();
	prev_curve_start_end.clear();
	prev_surf = "";
}

void SecondSkin::toolOnSetup(MEvent &)
{
	setHelpString(helpString);
	clear_quad_cache();
}

MStatus SecondSkin::doPress(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context)
{
	if (event.isModifierControl())
	{
		drawMode = EDDrawMode::kNormal;
	}
	else if (event.isModifierShift())
	{
		drawMode = EDDrawMode::kTangent;
	}
	else
	{
		drawMode = EDDrawMode::kDefault;
	}

	// Get the active 3D view.
	//
	view = M3dView::active3dView();
	update_anchors();

	//// Create an array to hold the lasso points. Assume no mem failures
	//maxSize = initialSize;
	//lasso = (coord*)malloc(sizeof(coord) * maxSize);


	coord start;
	event.getPosition(start.h, start.v);
	auto snap_anchor_index = do_snap(start.toMPoint());
	if (snap_anchor_index > 0 && snap_anchor_index < anchors.size())
	{
		first_anchor = anchors[snap_anchor_index];
		start.h = static_cast<short>(first_anchor->point_2D.x);
		start.v = static_cast<short>(first_anchor->point_2D.y);
		first_anchored = true;
	}
	
	
	//num_points = 1;
	//lasso[0] = min = max = start;
    stroke.clear();
    stroke.push_back(start);
    min = start;
    max = start;


	firstDraw = true;

	return MS::kSuccess;
}

MStatus SecondSkin::doDrag(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context)
// Add to the growing lasso
{

	if (!firstDraw) {
		//	Redraw the old lasso to clear it.
		draw_stroke(drawMgr);
	}
	else {
		firstDraw = false;
	}

	coord currentPos;
	event.getPosition(currentPos.h, currentPos.v);
	append_stroke(currentPos.h, currentPos.v);

	////	Draw the new lasso.
	draw_stroke(drawMgr);

	return MS::kSuccess;
}

MStatus SecondSkin::doRelease(MEvent & event, MHWRender::MUIDrawManager& drawMgr, const MHWRender::MFrameContext& context)
// Selects objects within the lasso
{
	MStatus							stat;
	MSelectionList					incomingList, boundingBoxList, newList;

	if (!firstDraw) {
		// Redraw the lasso to clear it.
		//view.beginXorDrawing(true, true, 1.0f, M3dView::kStippleDashed);
		draw_stroke(drawMgr);
	}
	//view.viewToObjectSpace

	bool first_point_known = false;
	bool last_point_known = false;
	MPoint first_world_point;
	MPoint last_world_point;

	if (drawing_quad)
	{
		if (prev_curves.size() >= 1)
		{
			first_point_known = true;
			first_world_point = prev_curve_start_end.back().second;
		}

		if (prev_curves.size() >= 3)
		{
			last_point_known = true;
			last_world_point = prev_curve_start_end.front().first;
		}
	}

	if (!last_point_known)
	{
		coord currentPos;
		event.getPosition(currentPos.h, currentPos.v);
		auto snap_anchor_index = do_snap(currentPos.toMPoint());
		if (snap_anchor_index > 0 && snap_anchor_index < anchors.size())
		{
			auto& acr = anchors[snap_anchor_index];
			append_stroke(static_cast<short>(acr->point_2D.x), static_cast<short>(acr->point_2D.y));
			last_point_known = true;
			last_world_point = acr->point_3D;
            last_anchor = acr;
		}
	}
	if (!first_point_known)
	{
		if (first_anchored)
		{
			first_point_known = true;
			first_world_point = first_anchor->point_3D;
		}
	}

	// get selection location
	MGlobal::getActiveSelectionList(incomingList);
	MItSelectionList iter(incomingList);
	//MPoint selectionPoint(0, 0, 0);

	MFnMesh * selected_mesh = nullptr;

	for (; !iter.isDone(); iter.next())
	{
		MDagPath dagPath;
		auto stat = iter.getDagPath(dagPath);

		if (stat)
		{

			if (dagPath.hasFn(MFn::kTransform))
			{
				MFnTransform fn(dagPath);
				//selectionPoint = fn.getTranslation(MSpace::kWorld);
				dagPath.extendToShape();
			}

			if (dagPath.hasFn(MFn::kMesh))
			{
				MStatus stat_mesh;
				selected_mesh = new MFnMesh(dagPath, &stat_mesh);

				if (stat_mesh)
				{
					break;
				}
				else
				{
					selected_mesh = nullptr;
					continue;
				}
			}
		}
	}

	// TODO: rebuild kd-tree only when view is changed
	rebuild_kd(selected_mesh);

	// generate curve
	auto tan_mode = drawMode == EDDrawMode::kTangent;
	auto norm_mode = drawMode == EDDrawMode::kNormal;
	MString new_curve = create_curve(stroke, selected_mesh, first_point_known, last_point_known, first_world_point, last_world_point, tan_mode, norm_mode);
	
	
	// generate surface or volume
	if (new_curve != "")
	{
		if (prev_curves.size() == 4) {
			std::string surface_command;
			surface_command.reserve(500);
			surface_command.append("proc string __ed_draw_surf() { \n");
			surface_command.append("string $slct[]=`ls- sl`;\n");
			surface_command.append("select -r ");
			std::list<std::string> curve_names;
			while (!prev_curves.empty())
			{
				curve_names.push_back(prev_curves.front().asChar());
				surface_command.append(prev_curves.front().asChar());
				prev_curves.pop_front();
				prev_curve_start_end.pop_front();
				surface_command.append(" ");
			}
			surface_command.append(";\n");
			surface_command.append("string $nurbssurf[] = `boundary -ch 1 -or 0 -ep 0 -rn 0 -po 0 -ept 0.01 ");
			while (!curve_names.empty())
			{
				surface_command.append(" ");
				surface_command.append("\"" + curve_names.front() + "\"");
				curve_names.pop_front();
			}
			surface_command.append("`;\n");
			surface_command.append("select $slct; \n");
			surface_command.append("return $nurbssurf[0]; \n } \n");
			surface_command.append("__ed_draw_surf();");
			MGlobal::executeCommand(MString(surface_command.c_str()), prev_surf);
			drawn_shapes.push_back(prev_surf);
		}

		if (norm_mode && prev_surf != "")
		{
			// extrude previous surface;
			float distance = (drawn_curves.back().start - drawn_curves.back().end).length();
			std::string extrude_command;
			extrude_command.reserve(500);

			extrude_command.append("string $mesh_surf[]= `nurbsToPoly -mnd 1  -ch 1 -f 1 -pt 1 -pc 200 -chr 0.9 -ft 0.01 -mel 0.001 -d 0.1 -ut 1 -un 3 -vt 1 -vn 3 -uch 0 -ucr 0 -cht 0.01 -es 0 -ntr 0 -mrt 0 -uss 1 \"");
			extrude_command.append(prev_surf.asChar());
			extrude_command.append("\"`;\n");
			extrude_command.append("polyExtrudeFacet -ltz " + std::to_string(distance) + " -constructionHistory 1 -keepFacesTogether 1 -divisions 4 -twist 0 -taper 1 -off 0 -thickness 0 -smoothingAngle 30 $mesh_surf[0];");
			MGlobal::executeCommand(MString(extrude_command.c_str()));
			//TODO: Move this out

			prev_curves.clear();
			prev_curve_start_end.clear();
		}
	}
	

	stroke.clear();
	first_anchored = false;
    first_anchor = nullptr;
    last_anchor = nullptr;

	if (selected_mesh)
		delete selected_mesh;
	selected_mesh = nullptr;

	return MS::kSuccess;
}

MStatus SecondSkin::drawFeedback(MHWRender::MUIDrawManager & drawMgr, const MHWRender::MFrameContext & context)
{
	//draw_stroke(drawMgr);
	return MStatus::kSuccess;
}

size_t SecondSkin::do_snap(const MPoint & input_end_point)
{
	const float radius = 5;

	float pt[] = { input_end_point.x , input_end_point.y, 0 };
	
	size_t out_index = 0;
	float out_dist_squared = 0;
	anchors_kd_2d->knnSearch(pt, 1, &out_index, &out_dist_squared);

	if (out_dist_squared < radius * radius && out_index > 0)
	{
		return out_index;
	}

	return -1;
}

void SecondSkin::completeAction()
{
}

void SecondSkin::deleteAction()
{
	if (drawn_shapes.size() == 0) return;
	MString delete_command;
	delete_command += "delete ";
	auto last = drawn_shapes.back();
	delete_command += last;
	drawn_shapes.pop_back();
	delete_command += ";";

	if (prev_curves.size() > 0 && prev_curves.back() == last)
	{
		// step back quad cache
		prev_curves.pop_back();
		prev_curve_start_end.pop_back();
	}

	if (drawn_curves.size() > 0 && drawn_curves.back().name == last)
	{
		drawn_curves.pop_back();
	}
	 
	MGlobal::executeCommand(delete_command);
}

MString SecondSkin::create_curve(std::vector<coord>& screen_points, MFnMesh* selected_mesh, bool start_known, bool end_known, const MPoint & start_point, MPoint & end_point, bool tangent_mode, bool normal_mode)
{
	if (!selected_mesh)
	{
		return MString();
	}

	auto num_points = screen_points.size();

	std::vector<MPoint> world_points;
	world_points.reserve(num_points);
	std::vector<bool> hit_list;
	hit_list.reserve(num_points);
	std::vector<std::pair<MPoint, MVector>> rays;
	rays.reserve(num_points);

	unsigned hit_count = 0;
	// calculate points in world space
	for (unsigned i = 0; i < num_points; i++)
	{
		MPoint ray_origin;
		MVector ray_direction;

		view.viewToWorld(screen_points[i].h, screen_points[i].v, ray_origin, ray_direction);
		rays.push_back(std::pair<MPoint, MVector>(ray_origin, ray_direction));
		MPoint world_point;
		bool hit = false;

			//MIntArray faceids;
			MFloatPoint hit_point;
			float hit_param;
			int hit_face;
			int hit_tri;
			float hit_bary1;
			float hit_bary2;

			bool intersected = selected_mesh->closestIntersection(ray_origin,
				ray_direction,
				nullptr,
				nullptr,
				true,
				MSpace::kWorld,
				10000, // maxParam
				false, // testBothDirections
				nullptr,
				hit_point,
				&hit_param,
				&hit_face,
				&hit_tri,
				&hit_bary1,
				&hit_bary2
				);

			if (intersected)
			{
				world_point = hit_point;
				hit = true;
				hit_count++;
			}
		

		world_points.push_back(world_point);
		hit_list.push_back(hit);
	}
	if (start_known)
	{
		world_points[0] = start_point;
	}
	if (end_known)
	{
		world_points[num_points - 1] = end_point;
	}

	if (world_points.size() > 2)
	{
		bool projecting_normal = false;
		if (hit_count == 0)
		{
			project_contour(screen_points, world_points, hit_list, selected_mesh, rays, start_known, end_known);
			// classify this curve as shell contour
			setHelpString("Classified: Shell Contour!");

			// TODO: SHAPE MATCHING!
		}
		else if ((is_normal(screen_points, world_points, hit_list, selected_mesh) || normal_mode) && (hit_list[0] || hit_list[num_points - 1]))
		{
			projecting_normal = true;
			project_normal(screen_points, world_points, hit_list, selected_mesh, rays, start_known, end_known);
			setHelpString("Classified: Normal!");
		}
		else if (tangent_mode)
		{
			// TODO: actually use is_tangent the same time as force tangent
			project_tangent(screen_points, world_points, hit_list, selected_mesh, rays, start_known, end_known);
			setHelpString("Classified: Tangent Plane!");
		}
		else
		{
			project_shell(screen_points, world_points, hit_list, selected_mesh, rays, start_known, end_known);
			setHelpString("Classified: Shell Projection!");
		}

		// create the curve
		std::string curve_command;
		curve_command.reserve(world_points.size() * 40 + 200);
		curve_command.append("proc string __ed_draw_curve() { \n");
		curve_command.append("string $slct[]=`ls- sl`;\n");
		curve_command.append("string $cv = `curve");
		for (auto & p : world_points)
		{
			curve_command.append(" -p ");
			curve_command.append(std::to_string(p.x));
			curve_command.append(" ");
			curve_command.append(std::to_string(p.y));
			curve_command.append(" ");
			curve_command.append(std::to_string(p.z));
		}
		curve_command.append("`;\n");
		// smooth the curve
		curve_command.append("rebuildCurve -ch 1 -rpo 1 -rt 0 -end 1 -kr 0 -kcp 0 -kep 1 -kt 0 -s 8 -d 3 -tol 0.01 $cv; \n");
		curve_command.append("select $slct; \n");
		curve_command.append("return $cv; \n } \n");
		curve_command.append("__ed_draw_curve();");

		MString curve_name;
		MGlobal::executeCommand(MString(curve_command.c_str()), curve_name);
		if (!projecting_normal)
		{
			prev_curves.push_back(curve_name);
			prev_curve_start_end.push_back(std::pair<MPoint, MPoint>(world_points[0], world_points[world_points.size() - 1]));
		}

        auto cv = DrawnCurve(world_points[0], world_points[world_points.size() - 1], curve_name);
        cv.start_anchor = first_anchor;
        cv.end_anchor = last_anchor;
        MPoint dummypoint2D;
        if (!cv.start_anchor)
        {
            cv.start_anchor.reset(new EDAnchor(dummypoint2D, cv.start));
        }
        if (!cv.end_anchor)
        {
            cv.end_anchor.reset(new EDAnchor(dummypoint2D, cv.end));
        }
		drawn_curves.push_back(cv);
        
		drawn_shapes.push_back(curve_name);
		return curve_name;

	}
	else
	{
		return MString();
	}
}



void SecondSkin::search_loop_from(DrawnCurve &cv, int current_depth, std::list<MString>& loop_list)
{
    for (DrawnCurve& cv2 : drawn_curves)
    {
        if (cv.name == cv2.name) continue;
        
        if (cv.end_anchor == cv2.start_anchor)
        {
            current_depth++;
            search_loop_from(cv2, current_depth, loop_list);
        }
    }
}

MString SecondSkin::create_surface_from_loop(DrawnCurve& cv)
{
	std::list<MString> loop;
	search_loop_from(cv, 0, loop);
	return MString();
}

bool SecondSkin::is_normal(const std::vector<coord> & screen_points, const std::vector<MPoint> & world_points, const std::vector<bool> & hit_list,
	    const MFnMesh * selected_mesh) const
{
	if (!selected_mesh)
	{
		return false;
	}

	if (hit_list[0])
	{
		// if starting point is normal
		//unsigned sample_num = 0;
		auto num_points = screen_points.size();

		MPoint p0 = screen_points[0].toMPoint();
		MVector current_tang;
		for (int i = 1; i < tang_samples; i++)
		{
			if (i >= num_points) break;
			current_tang += screen_points[i].toMPoint() - p0;
			//sample_num++;
		}
		current_tang.normalize();

		MPoint closest_point;
		MVector surface_normal;
		selected_mesh->getClosestPointAndNormal(world_points[0], closest_point, surface_normal, MSpace::kWorld);
		surface_normal.normalize();
		MPoint point_plus_normal = closest_point + surface_normal;

		coord cp0, cp1;
		view.worldToView(closest_point, cp0.h, cp0.v);
		view.worldToView(point_plus_normal, cp1.h, cp1.v);
		MVector norm_proj = cp1.toMPoint() - cp0.toMPoint();
		norm_proj.normalize();

		// * is dot... strange Maya...
        auto v = 1.0 - (current_tang * norm_proj);
		return v < normal_threshold;
	}

	// todo: last_hit

}

//bool SecondSkin::is_tangent()const
//{
//	// FIXME: curvature has problems.
//	// not very using this function since I am just using SHIFT to force tangent
//	auto num_points = stroke.size();
//	double menger_curvature = 0;
//    unsigned valid_points = 0;
//	for (int i = 0; i < num_points - 2; i++){
//        MPoint x = stroke[i].toMPoint();
//		MPoint y = stroke[i + 1].toMPoint();
//		MPoint z = stroke[i + 2].toMPoint();
//		MVector xy = x - y;
//		MVector zy = z - y;
//		MVector zx = z - x;
//		double area = (xy^zy).length();
//
//        if (xy.length() == 0 && zy.length() == 0 && zx.length() == 0){
//            continue;
//        }
//
//		if (xy.length() != 0 && zy.length() != 0 && zx.length() != 0){
//			menger_curvature += 4 * area / (xy.length()*zy.length()*zx.length());
//		}
//        valid_points += 1;
//	}
//
//    if (valid_points == 0) return false;
//
//	menger_curvature /= valid_points;
//	if (menger_curvature <= 0.2){
//		return true;
//	}
//	return false;
//}

void SecondSkin::project_normal(std::vector<coord> & screen_points, std::vector<MPoint>& world_points, const std::vector<bool>& hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known)
{
	if (!selected_mesh)
	{
		return;
	}

	if (hit_list[0])
	{

		MPoint closest_point;
		MVector surface_normal;
		selected_mesh->getClosestPointAndNormal(world_points[0], closest_point, surface_normal, MSpace::kWorld);
		surface_normal.normalize();

		auto normal = EDMath::minimumSkewViewplane(rays[0].second, surface_normal);
		auto point = world_points[0];
		auto point_num = rays.size();

		for (int i = 0; i < point_num; i++)
		{
			world_points[i] = EDMath::projectOnPlane(point, normal, rays[i].first, rays[i].second);
		}
	}
	// todo: last hit
}
///
// Find a point on a camera ray that is nearest to the mesh
///
MPoint SecondSkin::find_point_nearest_to_mesh(const MFnMesh * selected_mesh, const MPoint & ray_origin, const MVector & ray_direction, const coord & screen_coord, float & ret_height) const
{
	if (!selected_mesh)
	{
		return ray_origin;
	}

	float pt[] = { screen_coord.h , screen_coord.v, 0 };
	size_t out_index = 0;
	float out_dist_squared = 0;
	kd_2d->knnSearch(pt, 1, &out_index, &out_dist_squared);

	// nearest point (I am just using vertex for now) on the mesh
	auto temp = mesh_pts.pts[out_index];
	MPoint p_on_mesh(temp.x, temp.y, temp.z);

	auto dist = (ray_direction * (p_on_mesh - ray_origin));
	if (dist < 0)
	{
		return ray_origin;
	}

	auto p_on_ray = ray_origin + dist * ray_direction;

	ret_height = (p_on_ray - p_on_mesh).length();

	return p_on_ray;
}



void SecondSkin::project_contour(std::vector<coord> & screen_points, std::vector<MPoint>& world_points, const std::vector<bool>& hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>>& rays,bool first_point_known, bool last_point_known)
{
	if (!selected_mesh || !kd_2d || world_points.size() < 2)
	{
		return;
	}

	// TODO: with shape matching
	// TODO: find nearest point on mesh, not vertex

	auto length = rays.size();
	float dummy;
	auto s0 = find_point_nearest_to_mesh(selected_mesh, rays[0].first, rays[0].second, screen_points[0], dummy);
	auto sn = find_point_nearest_to_mesh(selected_mesh, rays[length - 1].first, rays[length - 1].second, screen_points[length - 1], dummy);

	if (first_point_known)
	{
		s0 = world_points[0];
	}
	if (last_point_known)
	{
		sn = world_points[length - 1];
	}

	if (s0.isEquivalent(sn)) return;

	auto d = (sn - s0).normal();
	auto normal = EDMath::minimumSkewViewplane(rays[0].second, d);

	world_points[0] = s0;
	world_points[length - 1] = sn;

	for (int i = 1; i < length - 1; i++)
	{
		world_points[i] = EDMath::projectOnPlane(s0, normal, rays[i].first, rays[i].second);
	}

}

double interpolate_height(const MPoint& p, const MPoint& p_start, const MPoint& p_end, double h_start, double h_end)
{
	auto w1 = (p - p_start).length();
	auto w2 = (p - p_end).length();

	return (w2 * h_start + w1 * h_end) / (w1 + w2);

}

///                 
// Shell Projection
///
void SecondSkin::project_shell(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known)
{
	// TODO: for intersection: find all points in distance 2h
	// TODO: get all triangles connected to the points by MFnMesh:getTriangles()
	// TODO: extrude them by h, connect them.
	// TODO: cast the ray again and find the intersection

	if (!selected_mesh || !kd_2d || world_points.size() < 2)
	{
		return;
	}

	// TODO: snaping to a known height and do interpolation

	auto length = rays.size();
	float start_height = 0, end_height = 0;
	//MPoint s0 = world_points[0];
	//MPoint sn = world_points[length - 1];
	//int iter_start = 0, iter_end = length;

	if (!hit_list[0] && !first_point_known)
	{
		world_points[0] = find_point_nearest_to_mesh(selected_mesh, rays[0].first, rays[0].second, screen_points[0], start_height);
	}
	start_height = static_cast<double>(EDMath::distance_to_mesh(selected_mesh, world_points[0]));

	if (!hit_list[length - 1] && !last_point_known)
	{
		world_points[length - 1] = find_point_nearest_to_mesh(selected_mesh, rays[length - 1].first, rays[length - 1].second, screen_points[length - 1], end_height);
	}

	end_height = static_cast<double>(EDMath::distance_to_mesh(selected_mesh, world_points[length - 1]));

	int first_miss = -1, last_miss = -1;
	for (size_t i = 1; i <= length - 2; i++)
	{
		if (!hit_list[i])
		{
			if (first_miss == -1)
				first_miss = i;
			last_miss = i;
		}
		else
		{
			// TODO: known heights other than h0 and hn.
			auto h = interpolate_height(rays[i].first, rays[0].first, rays[length - 1].first, start_height, end_height);
			world_points[i] = (-rays[i].second) * h + world_points[i];

			if (first_miss != -1 && last_miss != -1)
			{
                auto plane_normal = EDMath::minimumSkewViewplane(rays[first_miss - 1].second
                    , world_points[last_miss + 1] - world_points[first_miss - 1]);
				for (size_t j = first_miss; j <= last_miss; j++)
				{
                    world_points[j] = EDMath::projectOnPlane(world_points[first_miss - 1], plane_normal, rays[j].first, rays[j].second);
				}
			}
			first_miss = -1;
			last_miss = -1;
		}
	}
	if (first_miss != -1 && last_miss != -1)
	{
        auto plane_normal = EDMath::minimumSkewViewplane(rays[first_miss - 1].second
            , world_points[last_miss + 1] - world_points[first_miss - 1]);
        for (size_t j = first_miss; j <= last_miss; j++)
        {
            world_points[j] = EDMath::projectOnPlane(world_points[first_miss - 1], plane_normal, rays[j].first, rays[j].second);
        }
	}
}
// tangent projection
void SecondSkin::project_tangent(std::vector<coord> & screen_points, std::vector<MPoint> & world_points, const std::vector<bool> & hit_list, const MFnMesh * selected_mesh, std::vector<std::pair<MPoint, MVector>> & rays, bool first_point_known, bool last_point_known)
{
	if (!selected_mesh || !kd_2d || world_points.size() < 2) {
		return;
	}
	auto length = rays.size();

	//bool first_point_known = false, last_point_known = false;
	//if (drawing_quad)
	//{
	//	if (prev_curves.size() >= 1)
	//	{
	//		world_points[0] = prev_curve_start_end.back().second;
	//		//start_height = static_cast<double>(EDMath::distance_to_mesh(selected_mesh, world_points[0]));
	//		// TODO: height
	//		first_point_known = true;
	//	}

	//	if (prev_curves.size() >= 3)
	//	{
	//		world_points[length - 1] = prev_curve_start_end.front().first;
	//		//end_height = static_cast<double>(EDMath::distance_to_mesh(selected_mesh, world_points[length - 1]));
	//		last_point_known = true;
	//	}
	//}
	
	//determine the height of the tangent plane and the middle point on that plane
	//assume the average height is the height of the middle point
	float h = 0.0;
	int mid_index = int(length / 2);
	MPoint nearest_point = find_point_nearest_to_mesh(selected_mesh, rays[mid_index].first, rays[mid_index].second, screen_points[mid_index], h);
	MPoint middle_point = (-rays[mid_index].second)*h + world_points[mid_index];

	//project each stroke points on the base layer and get each normal
	MVector normal;
	MPoint closest_point;
	MVector sum_normal = MVector(0.0, 0.0, 0.0);

	for (int i = 0; i < length; i++) {
		selected_mesh->getClosestPointAndNormal(world_points[i], closest_point, normal, MSpace::kWorld);
		sum_normal += normal;
	}
	//calculate the average normal as the tangent plane normal
	sum_normal = MVector(sum_normal.x / length, sum_normal.y / length, sum_normal.z / length);
	MVector plane_normal = sum_normal.normal();
	//project all the point on to the tangent plane
	for (int i = 0; i < length; i++) {
		world_points[i] = (-rays[i].second) * h + world_points[i];
		world_points[i] = EDMath::projectOnPlane(middle_point, plane_normal, rays[i].first, rays[i].second);
	}
	
}

void SecondSkin::rebuild_kd(const MFnMesh * selected_mesh)
{
	mesh_pts.clear();
	mesh_pts_2d.clear();
	if (!selected_mesh)
	{
		kd_2d = nullptr;
		return;
	}
	MPointArray pts_array;
	selected_mesh->getPoints(pts_array, MSpace::kWorld);

	auto length = pts_array.length();

	mesh_pts.pts.resize(length);
	for (int i = 0; i < length; i++)
	{
		mesh_pts.pts[i].x = pts_array[i].x;
		mesh_pts.pts[i].y = pts_array[i].y;
		mesh_pts.pts[i].z = pts_array[i].z;
	}

	mesh_pts_2d.pts.resize(length);
	for (int i = 0; i < length; i++)
	{
		short x, y;
		view.worldToView(pts_array[i], x, y);
		mesh_pts_2d.pts[i].x = x;
		mesh_pts_2d.pts[i].y = y;
		mesh_pts_2d.pts[i].z = 0;
	}

	kd_2d.reset(new EDMath::KDTree2D(2 /*dim*/, mesh_pts_2d, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
	kd_2d->buildIndex();

}

void SecondSkin::append_stroke(short x, short y)
{
	// TODO: make this axis independent and can increment by length?
	int		cy, iy, ix, ydif, yinc, i;
	float	fx, xinc;

	auto num_points = stroke.size();

	iy = (int)stroke[num_points - 1].v;
	ix = (int)stroke[num_points - 1].h;
	ydif = abs(y - iy);
	if (ydif == 0)
		return;

	// Keep track of smallest rectangular area of the screen that
	// completely contains the stroke.
	if (min.h > x)
		min.h = x;
	if (max.h < x)
		max.h = x;
	if (min.v > y)
		min.v = y;
	if (max.v < y)
		max.v = y;

	if (((int)y - iy) < 0)
		yinc = -1;
	else
		yinc = 1;

	xinc = (float)((int)x - ix) / (float)ydif;
	fx = (float)ix + xinc;
	cy = iy + yinc;
	for (i = 0; i < ydif; i++) {
		coord new_coord;
		new_coord.h = (short)fx;
		new_coord.v = (short)cy;

		stroke.push_back(new_coord);

		fx += xinc;
		cy += yinc;
		num_points++;
	}

	return;
}

void SecondSkin::update_anchors()
{
	anchors.clear();
	anchors_2d.clear();
	for (auto & curve : drawn_curves)
	{
		short x, y;
		view.worldToView(curve.start, x, y);
        curve.start_anchor->point_2D = MPoint(x, y, 0);
		anchors.push_back(curve.start_anchor);
		anchors_2d.pts.push_back(EDMath::PointCloud<float>::Point(x, y, 0));
		view.worldToView(curve.end, x, y);
        curve.end_anchor->point_2D = MPoint(x, y, 0);
        anchors.push_back(curve.end_anchor);
		anchors_2d.pts.push_back(EDMath::PointCloud<float>::Point(x, y, 0));
	}
	anchors_kd_2d.reset(new EDMath::KDTree2D(2 /*dim*/, anchors_2d, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
	anchors_kd_2d->buildIndex();
}

void SecondSkin::draw_stroke(MHWRender::MUIDrawManager& drawMgr)
{
	auto num_points = stroke.size();
	MColor curve_color(0.1f, 0.13f, 0.95f);
	MColor anchor_color(0.9f, 0.1f, 0.5f);
	drawMgr.beginDrawable();
	drawMgr.setColor(curve_color);
	drawMgr.setLineWidth(3);
	for (unsigned i = 1; i < num_points; i++)
	{
		drawMgr.line2d(MPoint(stroke[i - 1].h, stroke[i - 1].v), MPoint(stroke[i].h, stroke[i].v));
	}
	drawMgr.setColor(anchor_color);
	drawMgr.setLineWidth(1);
	for (auto & anchor : anchors)
	{
		drawMgr.circle2d(anchor->point_2D, 4, false);
	}
	drawMgr.endDrawable();
}

//void SecondSkin::draw_anchors(MHWRender::MUIDrawManager & drawMgr)
//{
//}

MPoint coord::toMPoint() const
{
	return MPoint(static_cast<double>(h), static_cast<double>(v));
}
