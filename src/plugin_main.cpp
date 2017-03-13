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

#include <maya/MPxContext.h>
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <maya/MPxContextCommand.h>

#include "SecondSkin.h"


//////////////////////////////////////////////
// Command to create contexts
//////////////////////////////////////////////

class LassoContextCmd : public MPxContextCommand
{
public:
	LassoContextCmd() = default;
	virtual MPxContext* makeObj();
	static void*		creator();
};

MPxContext* LassoContextCmd::makeObj()
{
	return new SecondSkin;
}

void* LassoContextCmd::creator()
{
	return new LassoContextCmd;
}



//////////////////////////////////////////////
// plugin initialization
//////////////////////////////////////////////
MStatus initializePlugin(MObject obj)
{
	MStatus		status;
	MFnPlugin	plugin(obj, PLUGIN_COMPANY, "3.0", "Any");
	std::cout << "plugin loaded" << std::endl;
	status = plugin.registerContextCommand("lassoToolContext",
		LassoContextCmd::creator);

	if (!status) {
		status.perror("registerContextCommand");
		return status;
	}

	// Register User Interface
	MString cmd_create_ui = R"(source ")" + plugin.loadPath() + R"(/CreateUI.mel")";
	MString cmd_delete_ui = R"(source ")" + plugin.loadPath() + R"(/DeleteUI.mel")";
	MGlobal::executeCommand(cmd_create_ui);
	MGlobal::executeCommand(cmd_delete_ui);

	// set the mel scripts to be run when the plugin is loaded / unloaded
	status = plugin.registerUI("CreateUI", "DeleteUI");
	if (!status) {
		status.perror("registerUIScripts");
		return status;
	}

	return status;
}

MStatus uninitializePlugin(MObject obj)
{
	MStatus		status;
	MFnPlugin	plugin(obj);

	status = plugin.deregisterContextCommand("lassoToolContext");

	return status;
}
