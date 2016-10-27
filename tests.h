/*
 * tests.h
 *
 *  Created on: Oct 30, 2015
 *      Author: starside
 */

#ifndef TESTS_H_
#define TESTS_H_

#include "data.h"
#include "branchedChain.h"

#ifdef USE_GUI
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

vtkSmartPointer<vtkPolyData> pointsOnSphere(const int nLongitude, const int nAzimuth);
vtkSmartPointer<vtkPolyData> randomPointsOnSphere(const int numPoints);
vtkSmartPointer<vtkPolyData> showMonomers(const branchedChain *mol);

#endif

void testRng();
void testSaveState();
void testStructureFactor();
void testCenterOfMass();


#endif /* TESTS_H_ */
