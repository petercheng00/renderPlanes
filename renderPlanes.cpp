#include "utilities.h"
#include "RPEventHandler.h"
#include <osgDB/WriteFile>

string modelFile = "";
string mapFile = "";
string plyFile = "";
string iveFile = "";
string outputFile = "";

bool textureOn = true;
bool noTexture = false;
bool noSave = false;
bool showTriangles = false;
bool showConvex = false;

//liberally borrowed from modeling.exe
void doEarClipping( Geometry* planeGeometry, Vec2Array* planeVertices, bool makeConvex )
{
	vector<DrawElementsUInt*> triangles;
	int numPoints = planeVertices->size(); 

	vector<int> polyList; 
	for(int i = 0; i!=numPoints; ++i)
	{
		polyList.push_back(i); 
	}

	int curIndex = 0; 
	//iteratively find triangles with which to fill the polygon
	int curIterations = 0; 
	while(polyList.size() > 0)
	{
		//if there are only three vertices of the polygon left, we use them for the final triangle
		if(polyList.size() == 3)
		{
			DrawElementsUInt* lastTri = new DrawElementsUInt( PrimitiveSet::POLYGON, 0 );
			lastTri->push_back(polyList[0]);
			lastTri->push_back(polyList[1]);
			lastTri->push_back(polyList[2]);
			triangles.push_back(lastTri);
			polyList.clear();
		}
		//if there are more than three vertices left, we randomly choose a point and it's adjacent vertices to make a triangle
		else
		{	
			//if we've tried a lot of points, we just have to give up and make a triangle to allow the process to continue
			if(curIterations > numPoints + 10)
			{
				cerr << "Giving up after " << curIterations << " attempts " << polyList.size() << " remaining vertices" << endl;
				/*for (int k = 0; k < polyList.size(); k++) {
					cerr << polyList[k] << endl;
					cerr << planeVertices->at(polyList[k]).x() << ", " << planeVertices->at(polyList[k]).y() << endl;
				}*/
				DrawElementsUInt* currTriangle =
					new DrawElementsUInt( PrimitiveSet::POLYGON, 0 );
				currTriangle->push_back(polyList[0]);
				currTriangle->push_back(polyList[1]);
				currTriangle->push_back(polyList[2]);
				triangles.push_back(currTriangle);
				polyList.erase(polyList.begin() + 1);
				curIterations=0;
				continue;
			}


			//under normal circumstance, we choose a random point in the polygon and its adjacent vertices
			//int curIndex = rand()%polyList.size();
			curIndex = (curIndex + 1)%polyList.size();
			int prevIndex = curIndex - 1; 
			if(prevIndex < 0)
			{
				prevIndex = polyList.size()- 1; 
			}
			int nextIndex = curIndex + 1; 
			if(nextIndex > polyList.size() - 1)
			{
				nextIndex = 0; 
			}

			//for some lame reason there are a lot of duplicate points in some models.
			if (((planeVertices->at(polyList[prevIndex]).x() ==  planeVertices->at(polyList[curIndex]).x()) &&
				(planeVertices->at(polyList[prevIndex]).y() ==  planeVertices->at(polyList[curIndex]).y())) ||
				((planeVertices->at(polyList[nextIndex]).x() ==  planeVertices->at(polyList[curIndex]).x()) &&
				(planeVertices->at(polyList[nextIndex]).y() ==  planeVertices->at(polyList[curIndex]).y()))) {
				cerr << "removing duplicate points" << endl;
				DrawElementsUInt* currTriangle =
					new DrawElementsUInt( PrimitiveSet::POLYGON, 0 );
				currTriangle->push_back(polyList[prevIndex]);
				currTriangle->push_back(polyList[curIndex]);
				currTriangle->push_back(polyList[nextIndex]);
				triangles.push_back(currTriangle);
				polyList.erase(polyList.begin() + curIndex);
				curIterations=0;
				continue;
			}
				
				
			
			//now that we have a candidate triangle we have to perform a couple of checks

			//first, we need to make sure that no other points of the polygon are inside the current candidate triangle
			//bool anythingInside = false; 
			int numInside = 0;  // this is how we check if this plane was a ceiling or floor, we need a new test for walls
			for(int i = 0; i!=polyList.size(); ++i)
			{
				if(i == prevIndex || i == curIndex || i == nextIndex){continue;}
				Vec3 barycentric = calculateBaryCentricCoords(planeVertices->at(polyList[i]), planeVertices->at(polyList[prevIndex]), planeVertices->at(polyList[curIndex]), planeVertices->at(polyList[nextIndex])); 
				if(barycentric[0] > 0 && barycentric[0] < 1 && barycentric[1] > 0 && barycentric[1] < 1  && barycentric[2] > 0 && barycentric[2] < 1 )
				{
					numInside++; 
				}
			}			
			//second, we need to make sure the up direction of the candidate triangle matches with the up direction of the polygon
			Vec2 v1 = planeVertices->at(polyList[curIndex]) - planeVertices->at(polyList[prevIndex]); 
			Vec2 v2 = planeVertices->at(polyList[nextIndex]) - planeVertices->at(polyList[curIndex]); 
			Vec3 v1_3d = Vec3(v1[0], v1[1], 0);
			Vec3 v2_3d = Vec3(v2[0], v2[1], 0);
			double zCross = (v1_3d^v2_3d)[2];
			// we want zCross < 0 since we know points go clockwise

			//if both conditions are satisfied, we create the triangle and remove the current vertex from the polygon
			if(numInside == 0 && (zCross < 0))
			{
				DrawElementsUInt* currTriangle =
					new DrawElementsUInt( PrimitiveSet::POLYGON, 0 );
				currTriangle->push_back(polyList[prevIndex]);
				currTriangle->push_back(polyList[curIndex]);
				currTriangle->push_back(polyList[nextIndex]);
				triangles.push_back(currTriangle);
				polyList.erase(polyList.begin() + curIndex);
				curIterations=0;
				
			}
			else
			{
				//if the conditions are not met, we try again incrementing the curIteratiosn variable
				curIterations++;
			}
		}
	}
	//Done triangulating. Now combine triangles into convex polys. This is unnecessary for rendering.
	vector<DrawElementsUInt*> polygons = triangles;
	int currPolyInd = 0;
	bool update = true;
	int success = 0;
	while (update) {
		update = false;
		if (polygons.size() < 2){
			break;
		}
		for (int i = 0; i < polygons.size(); ++i){
			for (int j = i+1; j < polygons.size(); ++j){
				DrawElementsUInt* currPoly = polygons[i];
				DrawElementsUInt* otherPoly = polygons[j];
				DrawElementsUInt* newPoly = combinePolygons( currPoly, otherPoly);
				if (newPoly && isConvex(newPoly, planeVertices, 0)){
					success += 1;
					update = true;
					DrawElementsUInt* erase1 = polygons[i];
					DrawElementsUInt* erase2 = polygons[j];
					polygons.erase(polygons.begin() + i);
					if (i < j){
						polygons.erase(polygons.begin() + j - 1);
					}
					else{
						polygons.erase(polygons.begin() + j);
					}
					polygons.push_back(newPoly);
				}
			}
		}
	}
	cout << "Num pieces: " << polygons.size() << endl;
	for (int i = 0; i < polygons.size(); ++i){
		planeGeometry->addPrimitiveSet(polygons[i]);
	}
}




Vec2Array* convert3dTo2d( Vec3Array* vertices ) {
	Vec2Array* vertices_2d = new Vec2Array();
	for (int i = 0; i < vertices->size(); i++ ) {
		vertices_2d->push_back(Vec2(vertices->at(i)[0], vertices->at(i)[1]));
	}
	return vertices_2d;
}

void addLighting(Group* root) {
	ref_ptr<LightSource> lightSource = new LightSource;
	Light* light = lightSource->getLight();
	light->setLightNum(0);
	light->setPosition(osg::Vec4(0.0f,0.0f,0.0f,1.0f));
	light->setAmbient(osg::Vec4(0.7f, 0.7f,0.7f,1.0f));
	light->setDiffuse(osg::Vec4(0.0f, 0.0f,0.0f,1.0f));
	light->setSpecular(osg::Vec4(0.0f,0.0f,0.0f,1.0f));
	root->addChild(lightSource);
	StateSet* statelightON_OVRD = root->getOrCreateStateSet();
	statelightON_OVRD->setAttribute(light, StateAttribute::ON);
}

void applyColors(Group* root, bool showPieces){
	int numPlanes = root->getNumChildren();
	for (int i = 0; i < numPlanes; i++){
		osg::Vec4Array* colors = new osg::Vec4Array;
		Geode* currGeode = (Geode*)root->getChild(i);
		Geometry* currGeometry = (Geometry*)currGeode->getDrawable(0);
		int numPrimitives = currGeometry->getNumPrimitiveSets();
		for (int j = 0; j < numPrimitives; j++){
			int a = rand() % 1000;
			int b = rand() % 1000;
			int c = rand() % 1000;
			colors->push_back(osg::Vec4(float(a)/1000.0f,float(b)/1000.0f,float(c)/1000.0f, 1.0f));
			//float color = (i*numPrimitives+j) % (numPlanes * numPrimitives);
			//colors->push_back(osg::Vec4(color/(numPlanes*numPrimitives), color/(numPlanes*numPrimitives), color/(numPlanes*numPrimitives), 1.0f));
		}
		currGeometry->setColorArray(colors);
		if (showPieces) {
			currGeometry->setColorBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);
		}
		else {
			currGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
		}
	}
}

void applyTextures(Group* root, vector<string>& fileVect, vector<Vec2Array*>& coordVect) {

	for (int i = 0; i < fileVect.size(); i ++ ){
		if (fileVect[i] == ""){
			continue;
		}
		Geode* currGeode = (Geode*)root->getChild(i);
		Geometry* currGeometry = (Geometry*)currGeode->getDrawable(0);
		Vec2Array* currVectArray = coordVect[i];
		osg::Vec2Array* texcoords = new osg::Vec2Array(currVectArray->size());
		for (int j = 0; j < currVectArray->size(); j ++ ){
			(*texcoords)[j].set(currVectArray->at(j)[0], currVectArray->at(j)[1]);
		}
		currGeometry->setTexCoordArray(0,texcoords);

		osg::Texture2D* currTexture = new osg::Texture2D;
		osg::Image* currFace = osgDB::readImageFile(fileVect[i]);
		currTexture->setImage(currFace);
		currTexture->setDataVariance(osg::Object::DYNAMIC); 
		currTexture->setInternalFormatMode(osg::Texture::USE_IMAGE_DATA_FORMAT);
		currTexture->setWrap(osg::Texture2D::WRAP_T, osg::Texture::CLAMP_TO_EDGE); 
		currTexture->setWrap(osg::Texture2D::WRAP_S, osg::Texture::CLAMP_TO_EDGE); 
		currTexture->setWrap(osg::Texture2D::WRAP_R, osg::Texture::CLAMP_TO_EDGE); 
		currTexture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture::LINEAR); 
		currTexture->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture::LINEAR); 


		osg::StateSet* state = new osg::StateSet;
		state->setTextureAttributeAndModes(0, currTexture,osg::StateAttribute::ON);
		state->setMode(GL_CULL_FACE, osg::StateAttribute::OFF); 
		currGeode->setStateSet(state);
	}
}


void parseMapFile(string mapFile, vector<string>& fileVect, vector<Vec2Array*>& coordVect) {
	ifstream inFile(mapFile);
	if(!inFile.is_open()){
		cerr << "error opening .map file - no textures will be applied" << endl;
		return;
		}
	int numPlanes;
	inFile >> numPlanes;
	int i = 0;
	string currToken;
	double currX, currY;
	string imageFile;
	while (i < numPlanes) {
		inFile >> currToken;
		if (currToken == "SKIP_TO") {
			int target = numPlanes;
			inFile >> currToken;
			if (currToken != "END"){
				target = atoi(currToken.c_str());
			}
			while(i < target){
				fileVect.push_back("");
				coordVect.push_back(NULL);
				++i;
			}
		}
		else {
			int numCoords = atoi(currToken.c_str());
			inFile >> imageFile;
			fileVect.push_back(imageFile);
			Vec2Array* currCoords = new Vec2Array();
			for (int j = 0; j < numCoords; j++ ){
				inFile >> currX;
				inFile >> currY;
				currCoords->push_back(Vec2(currX, currY));
			}
			coordVect.push_back(currCoords);
			++i;
		}
	}


}

void parsePlyFile(Group* root)
{
	ifstream inFile(plyFile);
	if (!inFile.is_open()){ cerr << "error opening ply file" << endl;}
	
	int numVerts;
	int numTris;
	int numPlanes;
	Vec3Array* verts = new Vec3Array();

	string currLine;
	ifstream inStream(plyFile.c_str());
	getline(inStream, currLine);
	getline(inStream, currLine);
	//check ascii here, don't care for now
	//assuming all comments come here, might not be true
	while(currLine.find("end_header") != 0)
	{
		getline(inStream,currLine);
		if (currLine.find("element vertex") == 0)
		{
			numVerts = atoi(currLine.substr(14).c_str());
		}
		else if (currLine.find("element face") == 0)
		{
			numTris = atoi(currLine.substr(12).c_str());
		}
		else if (currLine.find("element region") == 0)
		{
			numPlanes = atoi(currLine.substr(14).c_str());
		}
	}
	int i1, i2, i3, i4;
	double v1, v2, v3;
	for (int i = 0; i < numVerts; ++i)
	{
		inStream >> v1 >> v2 >> v3;
		verts->push_back(Vec3(v1*1000, v2*1000, v3*1000));
	}

	for (int i = 0; i < numTris; ++i)
	{	
		inStream >> i1 >> i2 >> i3 >> i4;
		Geode* triGeode = new Geode();
		Geometry* triGeometry = new Geometry();
		triGeode->addDrawable( triGeometry );
		root->addChild( triGeode );
		Vec3Array* triVertices = new Vec3Array();
		triVertices->push_back( verts->at(i1) );
		triVertices->push_back( verts->at(i2) );
		triVertices->push_back( verts->at(i3) );
		triGeometry->setVertexArray(triVertices);
	}
}

void parseModelFile(Group* root, bool makeConvex) {
	ifstream inFile(modelFile); 
	if(!inFile.is_open()){cerr << "error opening file specifying planes" << endl;}
	int numPlanes;
	inFile >> numPlanes; 
	//iterate over each plane
	for(int i = 0; i!=numPlanes; ++i) {
		Geode* planeGeode = new Geode();
		Geometry* planeGeometry = new Geometry();
		planeGeode->addDrawable( planeGeometry );
		root->addChild( planeGeode );
		Vec3Array* planeVertices = new Vec3Array;
		//read in the number of delimiting points for a plane
		int numDelimitingPoints; 
		inFile >> numDelimitingPoints; 

		//read in the plane equation, which is comprised of the normal vector (nx,ny,nz) and d such that nx*x+ny*y+nz*z + d = 0
		double nx,ny,nz,d; 
		inFile >> nx >> ny >> nz >> d; 
		//if normal up or down, floor or ceiling, so set a mask
		if (nz == 1.0) {
			planeGeode->setNodeMask(0x00000001);
		}
		else if (nz == -1.0) {
			planeGeode->setNodeMask(0x00000002);
		}
		else {
			planeGeode->setNodeMask(0x00000004);
		}

		for(int j = 0; j!= numDelimitingPoints; ++j) {
			double vx,vy,vz; 
			inFile >> vx >> vy >> vz; 
			planeVertices->push_back( Vec3( 1000*vx, 1000*vy, 1000*vz) ); 
		}

		planeGeometry->setVertexArray( planeVertices );

		
		//dotransformation to make horizontal
		Vec3 normal = ((planeVertices->at(1) - planeVertices->at(0)) ^ (planeVertices->at(2) - planeVertices->at(1)));
		int v = 2;
		while (normal.length() == 0) {
			if (v == planeVertices->size()-1) {
				cerr << "Cannot calculate a normal vector for plane " << i << " stopping now" << endl;
				std::cin.get();
			}
			normal = ((planeVertices->at(v) - planeVertices->at(v-1)) ^ (planeVertices->at(v+1) - planeVertices->at(v)));
			v++;
		}
		Matrix rotateMat = Matrix();
		rotateMat.makeRotate(normal, Vec3(0, 0, 1));
		bool reverseFlip = false;
		for ( int j = 0; j < planeVertices->size() ; j++ ){
			planeVertices->at(j) = rotateMat * planeVertices->at(j);
		}
		float zVal = planeVertices->at(0)[2];


		Vec2Array* planeVertices2d = convert3dTo2d( planeVertices );

		//standardize so that points go clockwise
		//if sum > 0, then we are clockwise
		float sum = (planeVertices2d->at(0).x()-planeVertices2d->back().x())*(planeVertices2d->at(0).y()+planeVertices2d->back().y());
		for ( int j = 1; j < planeVertices2d->size(); j++) {
			sum += (planeVertices2d->at(j).x()-planeVertices2d->at(j-1).x())*(planeVertices2d->at(j).y()+planeVertices2d->at(j-1).y());
		}
		
		if (sum < 0) {
			//not sure how makerotate works, so doing it this safe way
			cerr << "flipping again so points go clockwise" << endl;
			Matrix unRotateMat = Matrix();
			unRotateMat.makeRotate(Vec3(0, 0, 1), normal);
			for ( int j = 0; j < planeVertices2d->size() ; j++ ){
				planeVertices->at(j)[0] = planeVertices2d->at(j)[0];
				planeVertices->at(j)[1] = planeVertices2d->at(j)[1];
				planeVertices->at(j)[2] = zVal;
				planeVertices->at(j) = unRotateMat * planeVertices->at(j);
			}

			rotateMat = Matrix();
			rotateMat.makeRotate(normal, Vec3(0, 0, -1));
			reverseFlip = true;
			for ( int j = 0; j < planeVertices->size() ; j++ ){
				planeVertices->at(j) = rotateMat * planeVertices->at(j);
			}
			zVal = planeVertices->at(0)[2];

			planeVertices2d = convert3dTo2d( planeVertices );
		}
		doEarClipping( planeGeometry, planeVertices2d, makeConvex );
		Matrix unRotateMat = Matrix();
		if (reverseFlip){
			unRotateMat.makeRotate(Vec3(0, 0, -1), normal);
		}
		else {
			unRotateMat.makeRotate(Vec3(0, 0, 1), normal);
		}
		for ( int j = 0; j < planeVertices2d->size() ; j++ ){
			planeVertices->at(j)[0] = planeVertices2d->at(j)[0];
			planeVertices->at(j)[1] = planeVertices2d->at(j)[1];
			planeVertices->at(j)[2] = zVal;
			planeVertices->at(j) = unRotateMat * planeVertices->at(j);
		}
	
	}
	inFile.close(); 
}


bool RPEventHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
    {
    switch(ea.getEventType())
    {
    case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            switch(ea.getKey())
            {			
			case 'f':
				camera->setCullMask(camera->getCullMask() ^ 0x00000001);
				return false;
				break;
			case 'c':
				camera->setCullMask(camera->getCullMask() ^ 0x00000002);
				return false;
				break;
            case 'w':
				camera->setCullMask(camera->getCullMask() ^ 0x00000004);
				return false;
				break;
			default:
				return false;
            } 
        }
    default:
        return false;
    }
}

int main(int argc, char** argv)
{
	if (argc > 1){
		ifstream inputFile(argv[1]);
		outputFile = string(argv[1]);
		outputFile = outputFile.substr(0, outputFile.find_last_of(".")) + ".ive";
		string fileName;
		getline(inputFile, fileName);
		while(inputFile){
			size_t l = fileName.find_last_of("."); 
			string extension = fileName.substr(l+1);
			if(extension == "model") {
				modelFile = fileName;
				cout << "Reading model file:" << endl;
				cout << fileName << endl;
			}
			else if(extension == "ply") {
				plyFile = fileName;
				cout << "Reading ply file:" << endl;
				cout << fileName << endl;
			}
			else if (extension == "map") {
				mapFile = fileName;
				cout << "Reading map file:" << endl;
				cout << fileName << endl;
			}
			else if (extension == "ive") {
				iveFile = fileName;
				cout << "Adding external ive file:" << endl;
				cout << fileName << endl;
			}
			else if (fileName == "noTexture") {
				noTexture = true;
			}
			else if (fileName == "noSave") {
				noSave = true;
			}
			else if (fileName == "showTriangles") {
				showTriangles = true;
			}
			else if (fileName == "showConvex") {
				showConvex = true;
			}
			getline(inputFile, fileName);
		}
	}
	else {
		cout << "No input files, quitting" << endl;
		return 0;
	}
    osgViewer::Viewer viewer;
    Group* root = new Group();
    vector<string> planeToImageFile;
    vector<Vec2Array*> planeToImageCoords;
	parseMapFile(mapFile, planeToImageFile, planeToImageCoords);
	if (plyFile == "")
	{
		parseModelFile(root, showConvex);
	}
	else
	{
		parsePlyFile(root);
	}
	if (noTexture){
		applyColors(root, showTriangles || showConvex);
	}
	else{
	    applyTextures(root, planeToImageFile, planeToImageCoords);
	}
	addLighting(root);
	if (iveFile != ""){
		osg::ref_ptr<Node> otherModel = osgDB::readNodeFile(iveFile);
		root->addChild(otherModel);
	}

	if (!noSave){
		cout << "Saving ive file:" << endl;
		cout << outputFile << endl;
		osgDB::writeNodeFile(*root, outputFile.c_str());
	}
	else {

		RPEventHandler * myRPEventHandler = new RPEventHandler(viewer.getCamera());
		viewer.addEventHandler(myRPEventHandler);

	 	viewer.setUpViewOnSingleScreen(0);
		viewer.setSceneData( root );
		viewer.setCameraManipulator(new osgGA::OrbitManipulator());
		viewer.realize();
 
        while( !viewer.done() )
        {
            viewer.frame();
        }
	}
    return 0;
}
