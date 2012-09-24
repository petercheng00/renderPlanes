#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <vector>

#include <osg/Node>
#include <osg/Group>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Texture2D>
#include <osgDB/ReadFile> 
#include <osgViewer/Viewer>
#include <osg/PositionAttitudeTransform>
#include <osgGA/TrackballManipulator>


using namespace std;
using namespace osg;


bool isConvex( DrawElementsUInt* poly, Vec2Array* vertices, int concave = 0);
DrawElementsUInt* combinePolygons( DrawElementsUInt* poly1, DrawElementsUInt* poly2);
bool sameSide( Vec2 p1, Vec2 p2, Vec2 a, Vec2 b );
bool vertexInsidePolygon( Vec2 testPoint, Vec2Array* vertices );
bool verticesInTriangle( int a, int b, int c, Vec2Array* vertices, vector<int>& valid );
Vec3 calculateBaryCentricCoords(Vec2& in, Vec2& v1, Vec2& v2, Vec2& v3);

#endif