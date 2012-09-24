#include "utilities.h"

bool isConvex( DrawElementsUInt* poly, Vec2Array* vertices, int concave ){
	//int threshhold = 5000.0;
	//int threshhold = 0;
	bool positive = false;
	for (int i = 0; i < poly->size(); i++){
		int j = (i+1) % poly->size();
		int k = (i+2) % poly->size();
		Vec2 e1 = vertices->at(poly->at(j)) - vertices->at(poly->at(i));
		Vec2 e2 = vertices->at(poly->at(k)) - vertices->at(poly->at(j));
		Vec3 e1_3d = Vec3(e1[0], e1[1], 0);
		Vec3 e2_3d = Vec3(e2[0], e2[1], 0);
		Vec3 cross = e1_3d ^ e2_3d;
		// do angles instead
		if (i == 0){
			positive = cross[2] > 0;
		}
		if ((cross[2] > 0) && !positive){
			float dist = (vertices->at(poly->at(k)) - vertices->at(poly->at(i))).length();
			//if (dist > threshhold){
				concave -= 1;
			//}
		}
		if ((cross[2] < 0) && positive){
			float dist = (vertices->at(poly->at(k)) - vertices->at(poly->at(i))).length();
			//if (dist > threshhold){
				concave -= 1;
			//}
		}
		if (concave < 0){
			return false;
		}
	}
	return true;
}

DrawElementsUInt* combinePolygons( DrawElementsUInt* poly1, DrawElementsUInt* poly2){
	DrawElementsUInt* candidatePoly = new DrawElementsUInt( PrimitiveSet::POLYGON, 0 );
	for (int i = 0; i < poly1->size(); i++){
		candidatePoly->push_back(poly1->at(i));
	}
	//Look for a shared edge
	int sharedIndex1 = -1;
	int sharedIndex2 = -1;
	bool opposite = false;
	for (int i = 0; i < candidatePoly->size(); i++) {
		if (sharedIndex1 != -1) {
			break;
		}
		for (int j = 0; j < poly2->size(); j++) {
			int ind11 = candidatePoly->at(i);
			int ind12 = candidatePoly->at((i+1) % candidatePoly->size());
			int ind21 = poly2->at(j);
			int ind22 = poly2->at((j+1) % poly2->size());
			if ((ind11 == ind21) && (ind12 == ind22)){
				sharedIndex1 = i;
				sharedIndex2 = j;
				break;
			}
			if ((ind11 == ind22) && (ind12 == ind21)){
				sharedIndex1 = i;
				sharedIndex2 = j;
				opposite = true;
				break;
			}
		}
	}
	if (sharedIndex1 == -1){
		return NULL;
	}
	//insert places elements before selected location, so we go backwards when ordering is the same.
	if (!opposite){
		for (int i = poly2->size()-1; i >= 0; --i) {
			int index = (sharedIndex2 + i + 2) % poly2->size();
			candidatePoly->insert(candidatePoly->begin() + sharedIndex1 + 1, poly2->at(index));
		}
	}
	else {
		int offset = sharedIndex1 + 1;
		for (int i = 0; i < poly2->size() - 2; ++i) {
			int index = (sharedIndex2 + i + 2) % poly2->size();
			candidatePoly->insert(candidatePoly->begin() + offset, poly2->at(index));
			offset += 1;
		}
	}
	return candidatePoly;

}

bool sameSide( Vec2 p1, Vec2 p2, Vec2 a, Vec2 b ) {
	Vec3 p1_3d = Vec3(p1[0], p1[1], 0);
	Vec3 p2_3d = Vec3(p2[0], p2[1], 0);
	Vec3 a_3d = Vec3(a[0], a[1], 0);
	Vec3 b_3d = Vec3(b[0], b[1], 0);

	Vec3 cp1 = ( b_3d - a_3d ) ^ ( p1_3d - a_3d );
	Vec3 cp2 = ( b_3d - a_3d ) ^ ( p2_3d - a_3d );
	return ( cp1 * cp2 >= 0 );
}

bool vertexInsidePolygon( Vec2 testPoint, Vec2Array* vertices ){
	int counter = 0;
	Vec2 p1 = vertices->front();
	for ( int i = 1; i <= vertices->size(); i++ ) {
		Vec2 p2 = vertices->at( i % vertices->size() );
		if ( testPoint[1] > min(p1[1], p2[1]) &&
			 testPoint[1] <= max(p1[1], p2[1]) &&
			 testPoint[0] <= max(p1[0], p2[0]) &&
			 p1[1] != p2[1]) {
				float xintercept = (testPoint[1]-p1[1])*(p2[0]-p1[0])/(p2[1]-p1[1])+p1[0];
				if (p1[0] == p2[0] || testPoint[0] <= xintercept) {
					counter++;
				}
		}
		p1 = p2;
	}
	return (counter % 2 == 1);
}

bool verticesInTriangle ( int a, int b, int c, Vec2Array* vertices, vector<int>& valid ) {
	for ( int i = 0; i < valid.size(); i++ ) {
		if ( valid[i] == a || valid[i] == b || valid[i] == c ) continue;
		if ( sameSide(vertices->at( valid[i] ), vertices->at(a), vertices->at(b), vertices->at(c)) &&
		     sameSide(vertices->at( valid[i] ), vertices->at(b), vertices->at(a), vertices->at(c)) &&
			 sameSide(vertices->at( valid[i] ), vertices->at(c), vertices->at(a), vertices->at(b)) ) {
				 return true;
		}
	}
	return false;
}

Vec3 calculateBaryCentricCoords(Vec2& in, Vec2& v1, Vec2& v2, Vec2& v3)
{
	//vec3 inPos = in->getPosition(); 
	double x = in[0]; 
	double y = in[1]; 

	//vec3 v1Pos = v1->getPosition(); 
	double xa = v1[0]; 
	double ya = v1[1]; 

	//vec3 v2Pos = v2->getPosition(); 
	double xb = v2[0]; 
	double yb = v2[1]; 

	//vec3 v3Pos = v3->getPosition(); 
	double xc = v3[0]; 
	double yc = v3[1]; 

	double gamma = ((ya-yb)*x + (xb-xa)*y + xa*yb - xb*ya)  /  ((ya-yb)*xc +(xb-xa)*yc + xa*yb - xb*ya);
	double beta = ((ya-yc)*x + (xc-xa)*y + xa*yc - xc*ya)  /  ((ya-yc)*xb + (xc-xa)*yb + xa*yc - xc*ya);
	double alpha = 1- beta - gamma;

	return Vec3(gamma, beta, alpha); 
}