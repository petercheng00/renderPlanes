#ifndef RPEVENTHANDLER_H
#define RPEVENTHANDLER_H
#include "utilities.h"

class RPEventHandler: public osgGA::GUIEventHandler
{
    public:
		RPEventHandler(Camera* c) : osgGA::GUIEventHandler() {
			camera = c;
		}
        virtual bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter&);
        virtual void accept(osgGA::GUIEventHandlerVisitor& v)   { v.visit(*this); };
		Camera* camera;
};

#endif