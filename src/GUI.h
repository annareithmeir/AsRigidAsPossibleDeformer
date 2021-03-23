#pragma once

#include <iostream>
#include <cstring>
#include <fstream>

#include "Eigen.h"
#include <unsupported/Eigen/src/MatrixFunctions/MatrixSquareRoot.h>
#include "SimpleMesh.h"
#include "ArapDeformer.h"

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/unproject.h>

#include <list>

class GUI {
public:
	//handles all rendering, camera movement, mesh interaction
	GUI(std::string filenameMesh, int iter, int get_weight_type, int get_estimation_type) {
		weight_type = get_weight_type;
		estimation_type = get_estimation_type;
		num_iterations = iter;
		displayMesh(filenameMesh);
	}

	std::set<int> fixedFaces, fixedFacesPreviousInit;
	std::set<int> handles, handlesPreviousInit;
	std::set<int> fixedVertices, fixedVerticesPreviousInit;
	int num_iterations;
	int weight_type;
	int estimation_type;

private:

	bool arap_mode = false;
	bool handle_selection_mode = false;
	bool faceSelection = false;
	bool mouse_down = false;
	bool handle_down = false;
	bool vertex_hit = false;
	bool arap_initialized = false;
	bool arap_running = false;
	bool deformer_initiated = false;

	int current_mouse_button = 0;
	int current_moving_handle = -1;
	int prev_moving_handle = -1;

	SimpleMesh source_mesh;
	ArapDeformer deformer;
	std::string meshName;

	Eigen::MatrixXd vertices, colors, handle_rep;
	Eigen::MatrixXi faces;

	//inspired by https://github.com/libigl/libigl/blob/master/tutorial/708_Picking/main.cpp

	void displayMesh(std::string filenameMesh) {
		//load mesh
		igl::readOFF(filenameMesh, vertices, faces);
		meshName = filenameMesh;

		//initialize our mesh structures based on loaded mesh
		if (!source_mesh.loadMesh(filenameMesh)) {
			std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
			return;
		}

		//init white face colors
		colors = Eigen::MatrixXd::Constant(faces.rows(), 3, 1);

		bool arapModeRef = arap_mode;
		bool handleSelectionModeRef = handle_selection_mode;

		igl::opengl::glfw::Viewer viewer;
		
		bool paint = true;
		//Called when mouse is pressed
		viewer.callback_mouse_down = [this, &paint](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (!arap_mode) {
				current_mouse_button = button;
				mouse_down = true;

				return selectionHandler(viewer, button);
			}
			else {
				current_mouse_button = button;
				mouse_down = true;

				int fid;
				Eigen::Vector3f bc;

				double x = viewer.current_mouse_x;
				double y = viewer.core().viewport(3) - viewer.current_mouse_y;

				if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
					viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
				{
					vertex_hit = true;
					int handleId = getClosestVertexIdFromBC(fid, bc);
					if (handles.find(handleId) != handles.end()) {
						handle_down = true;
						current_moving_handle = handleId;
						viewer.data().clear_points();
						set<int>::iterator itr;
						for (itr = fixedVertices.begin(); itr != fixedVertices.end(); itr++) {
							viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 0, 1));
						}
						requestArapInit();
						return true;
					}
				}
			}
			return false;
		};

		//Called when mouse is released
		viewer.callback_mouse_up = [this](igl::opengl::glfw::Viewer& viewer, int button, int) -> bool
		{
			if (button == current_mouse_button) {
				if (arap_mode && handle_down) {
					set<int>::iterator itr;
					for (itr = handles.begin(); itr != handles.end(); itr++) {
						viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
					}
				}
				mouse_down = false;
				vertex_hit = false;
				prev_moving_handle = current_moving_handle;
				handle_down = false;
			}
			return false;
		};

		//Called while the mouse is moving
		viewer.callback_mouse_move = [this](igl::opengl::glfw::Viewer& viewer, int, int) -> bool
		{
			//checks if mouse key is pressed + we are not in arap mode + have selected a vertex with the mouse key press
			if (mouse_down && vertex_hit) {

				if (!arap_mode) {
					return selectionHandler(viewer, current_mouse_button);
				}
				else {
					if (handle_down && !arap_running) {
						return displacementHandler(viewer);
					}
				}
			
			}
			return false;
		};

		//Called each time a key is pressed.
		viewer.callback_key_down = [this](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
		{
			//Fixed Vertex selection
			if (key == '1') {
				faceSelection = false;
				handle_selection_mode = false;
				arap_mode = false;
				for each (int fid in fixedFaces)
				{
					updateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				repaintVertices(viewer);
				std::string out = "Fixed Vertices Selection active, ARAP Paused";
				std::cout << out << std::endl;
				return true;
			}
			if (key == '2') {
				faceSelection = true;
				handle_selection_mode = false;
				arap_mode = false;
				for each (int fid in fixedFaces)
				{
					updateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				repaintVertices(viewer);
				std::string out = "Fixed Faces Selection active, ARAP Paused";
			    std::cout << out << std::endl;
				return true;
			}
			//Toggle handle selection
			if (key == '3') {
				faceSelection = false;
				arap_mode = false;
				for each (int fid in fixedFaces)
				{
					updateColor(fid, Eigen::Vector3d(1, 0, 0), viewer);
				}
				repaintVertices(viewer);
				handle_selection_mode = true;		
				std::cout << "Handle Selection Active, ARAP Paused" << std::endl;
			}
			//ARAP Mode
			if (key == '4') {
				faceSelection = false;
				handle_selection_mode = false;
				arap_mode = true;
				for each (int fid in fixedFaces)
				{
					updateColor(fid, Eigen::Vector3d(0, 0, 1), viewer);
				}
				//Repaint Vertices in blue
				viewer.data().clear_points();
				set<int>::iterator itr;
				for (itr = handles.begin(); itr != handles.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
				}
				for (itr = fixedVertices.begin(); itr != fixedVertices.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 0, 1));
				}
				std::string out = "ARAP active";
				std::cout << out << std::endl;
				return true;
			}
			if (key == '5') {
				viewer.open_dialog_save_mesh();
				return true;
			}
			return true;
		};

		//display mesh
		viewer.data().set_mesh(vertices, faces);
		viewer.data().point_size = 10;
		viewer.data().set_colors(colors);
		viewer.data().show_lines = false;
		viewer.data().double_sided = true;
		
		viewer.launch();
	}

	//Sends desired handle position to ARAP and updates vertices' positions
	void performARAP(Eigen::Vector3f handlePos, igl::opengl::glfw::Viewer& viewer) {
			arap_running = true;
			deformer.applyDeformation(getFixedVerticesFromFaces(), current_moving_handle, handlePos.cast<double>(), num_iterations); //flip flop optimization
			std::vector<Vertex> deformedVertices = deformer.m_mesh.getDeformedVertices();
			source_mesh = deformer.m_mesh;

			for (int i = 0; i < vertices.rows(); i++) {
				vertices.row(i) = deformedVertices[i].position;
			}

			//recompute normals for lighting
			viewer.data().compute_normals();
			//redraw, there appears to be a bottleneck in this method
			viewer.data().set_mesh(vertices, faces);

			arap_running = false;
	}

	//check if fixed or handle vertices have been added or removed or the active handle is another one than before -> re-init structs
	void requestArapInit() {
		if (handles.size() != 0 && (fixedFaces != fixedFacesPreviousInit || fixedVertices != fixedVerticesPreviousInit || handles != handlesPreviousInit || prev_moving_handle != current_moving_handle)) {
			arap_initialized = false;
			source_mesh = SimpleMesh();
			
			//reload potentially changed mesh from last process
			if (!source_mesh.loadMeshFromGUI(vertices, faces, getFixedVerticesFromFaces())) {
				std::cout << "Mesh file wasn't read successfully at location: " << meshName << std::endl;
				return;
			}

			//reinit to update weights + system matrix for current deformed mesh after previous completed interaction 
			deformer = ArapDeformer(&source_mesh, weight_type, estimation_type);
			deformer_initiated = true;
		
			fixedFacesPreviousInit = fixedFaces;
			handlesPreviousInit = handles;
			fixedVerticesPreviousInit = fixedVertices;
			arap_initialized = true;
			arap_running = false;
		}
	}

	//Retrieves all ids of fixed points from the chosen faces, includes the handle without duplicates
	std::vector<int> getFixedVerticesFromFaces() {
		std::set<int> fixedVerticesLocal;
		fixedVerticesLocal.clear();
		fixedVerticesLocal = fixedVertices;

		for (int face : fixedFaces) {
			for (int i = 0; i < 3; i++) {
				fixedVerticesLocal.insert(faces.row(face)(i));
			}
		}
		//if the fixed faces do not contain the prev handle -> make prev handle non-fixed!
		if (fixedVerticesLocal.find(current_moving_handle) != fixedVerticesLocal.end()) fixedVerticesLocal.erase(current_moving_handle);

		std::vector<int> fixedVerticesAsVector(fixedVerticesLocal.size());
		std::copy(fixedVerticesLocal.begin(), fixedVerticesLocal.end(), fixedVerticesAsVector.begin());
		fixedVerticesAsVector.push_back(current_moving_handle);
		return fixedVerticesAsVector;
	}

	// Calculates Handle Position in World Space according to mouse position
	bool displacementHandler(igl::opengl::glfw::Viewer& viewer) {
		//only run algorithm if initialized and not already running
		if (arap_initialized && !arap_running) {
			double x = viewer.current_mouse_x;
			double y = viewer.core().viewport(3) - viewer.current_mouse_y;

			Eigen::Vector3f handlePos = vertices.row(current_moving_handle).cast<float>();
			//Convert depth to view
			Eigen::Vector3d projection = igl::project(handlePos, viewer.core().view, viewer.core().proj, viewer.core().viewport).cast<double>();

			//Convert mouse position into world position
			Eigen::Vector3f worldPos = igl::unproject(Eigen::Vector3f(x, y, (float)projection.z()), viewer.core().view, viewer.core().proj, viewer.core().viewport);

			//send intended handle position to ARAP
			performARAP(worldPos, viewer);

			return true;
		}
		return false;
	}

	// Selects or Deselect faces/vertices on mouse click
	bool selectionHandler(igl::opengl::glfw::Viewer& viewer, int& mouseID) {
		int fid;
		Eigen::Vector3f bc;
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;
		double y = viewer.core().viewport(3) - viewer.current_mouse_y;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
			viewer.core().proj, viewer.core().viewport, vertices, faces, fid, bc))
		{
			vertex_hit = true;
			if(!faceSelection) fid = getClosestVertexIdFromBC(fid, bc);
			if (mouseID == 0) {
				if (handle_selection_mode) {
					select(fid, bc, handles, viewer);
				}
				else {
					if (faceSelection) {
						select(fid, bc, fixedFaces, viewer);
					}
					else {
						select(fid, bc, fixedVertices, viewer);
					}
				}
			}
			else if(mouseID == 2){
				if (handle_selection_mode) {
					unselect(fid, bc, handles, viewer);
				}
				else {
					if (faceSelection) {
						unselect(fid, bc, fixedFaces, viewer);
					}
					else {
						unselect(fid, bc, fixedVertices, viewer);
					}
				}
			}
			else {
				return false;
			}
			return true;
		}
		return false;
	}

	//Unselection utility function
	void unselect(int fid, Eigen::Vector3f& bc, std::set<int>& toSelect, igl::opengl::glfw::Viewer& viewer) {

		if (toSelect.find(fid) != toSelect.end()) {
			toSelect.erase(fid);
			//repaint
			if (!faceSelection) {
				viewer.data().clear_points();
				set<int>::iterator itr;
				for (itr = handles.begin(); itr != handles.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
				}
				for (itr = fixedVertices.begin(); itr != fixedVertices.end(); itr++) {
					viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(1, 0, 0));
				}
			}
			else {
				updateColor(fid, Eigen::Vector3d(1, 1, 1), viewer);
			}
		}
	}

	//Selection utility function
	void select(int fid, Eigen::Vector3f& bc, std::set<int>& toSelect, igl::opengl::glfw::Viewer& viewer)
	{
		auto newColor = handle_selection_mode ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);

		if (!faceSelection) {
			if (handles.find(fid) == handles.end() && fixedVertices.find(fid) == fixedVertices.end()) {
				int selectedVertex = fid;
				toSelect.insert(selectedVertex);
				viewer.data().add_points(vertices.row(selectedVertex), newColor.transpose());
				//checks that vertex is hit
				vertex_hit = true;
			}
		}
		else {
			//fixed faces are defined
			if (toSelect.find(fid) == toSelect.end()) {
				//select face
				toSelect.insert(fid);
				updateColor(fid, newColor, viewer);
			}
		}
	}

	//Updates overlay positions
	void repaintVertices(igl::opengl::glfw::Viewer& viewer) {
		viewer.data().clear_points();
		set<int>::iterator itr;
		for (itr = handles.begin(); itr != handles.end(); itr++) {
			viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(0, 1, 0));
		}
		for (itr = fixedVertices.begin(); itr != fixedVertices.end(); itr++) {
			viewer.data().add_points(vertices.row(*itr), Eigen::RowVector3d(1, 0, 0));
		}
	}

	//Finds the triangle vertex that is closest to a click on a face
	int getClosestVertexIdFromBC(int fid, Eigen::Vector3f& bc) {
		
		Eigen::Vector3i triangleVertices = faces.row(fid);		
		Eigen::MatrixXd closestPoints (3,3);

		for (int i = 0; i < 3; i++) {
			closestPoints.row(i) = vertices.row(triangleVertices(i));
		}
	
		Eigen::Vector3d sqdistances;
		//retrieve raycast result coordinates from barycentric ones
		Eigen::Vector3d queryPoint = bc(0) * closestPoints.row(0) + bc(1) * closestPoints.row(1)+ bc(2) * closestPoints.row(2);
		
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3d diff = closestPoints.row(i) - queryPoint.transpose();
			sqdistances(i) = std::sqrt(diff.dot(diff));
		}

		float min = 100;
		int minId = 0;
		float current;
		for (int i = 0; i < 3; i++) {
			current = sqdistances(i);
			if (min > current) {
				min = current;
				minId = i;
			}
		}

		return triangleVertices(minId);
		
	}

	void updateColor(int faceID, Eigen::Vector3d newColor, igl::opengl::glfw::Viewer& viewer) {
		colors.row(faceID) << newColor(0), newColor(1), newColor(2);
		viewer.data().set_colors(colors);
	}
};

//TODO make a nice menu with dropdown or buttons for modes, possibly brush size etc
class Menu {
public:
	Menu(igl::opengl::glfw::Viewer viewer)
	{
		this->viewer = viewer;
	}
private:
	igl::opengl::glfw::Viewer viewer;
	
};
