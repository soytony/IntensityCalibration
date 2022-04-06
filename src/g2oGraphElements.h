// This software is licenced under the Creative Commons
// (Attribution-NonCommercial-ShareAlike):
// http://creativecommons.org/licenses/by-nc-sa/4.0

// Author: Bastian Steder <steder@informatik.uni-freiburg.de>


#ifndef REMISSION_LEARNER_G2O_GRAPH_ELEMENTS_H
#define REMISSION_LEARNER_G2O_GRAPH_ELEMENTS_H

#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/core/base_binary_edge.h"
#include "g2o/core/base_multi_edge.h"
#include "g2o/core/hyper_graph_action.h"
//#include "aislib/stuff/macros.h"

//#include <QApplication>
//#include "ais3dTools/visualization/aluGLViewer/alu_glviewer.h"
#include "remissionCalibrationHelper.h"
#include "remissionCalibrationResult.h"

class RemissionLearnerVertex : public g2o::BaseVertex<1, double>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef BaseVertex<1, double> BaseClass;
    
    RemissionLearnerVertex();

    virtual void setToOriginImpl() {
      _estimate = 1.0;
    }

    virtual void oplusImpl(const double* update)
    {
      _estimate += *update;
    }
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
};

class RemissionLearnerUnaryEdge : public g2o::BaseUnaryEdge<1, double, RemissionLearnerVertex>
{
  public:
    typedef g2o::BaseUnaryEdge<1, double, RemissionLearnerVertex> BaseClass;
    
    RemissionLearnerUnaryEdge();
    
    static inline double jacobian(double vertexValue);
    inline void computeError();
    inline void linearizeOplusNumeric();
    inline void linearizeOplusJacobian();
    virtual void linearizeOplus();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class RemissionLearnerPriorEdge : public g2o::BaseUnaryEdge<1, double, RemissionLearnerVertex>
{
  public:
    typedef g2o::BaseUnaryEdge<1, double, RemissionLearnerVertex> BaseClass;
    
    RemissionLearnerPriorEdge(float value=1.0f, float weight=1.0f);
    
    static inline double jacobian(double vertexValue);
    inline void computeError();
    virtual void linearizeOplus();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};


class RemissionLearnerMultiEdge_n : public g2o::BaseMultiEdge<1, double> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    typedef g2o::BaseMultiEdge<1, double> BaseClass;
    
    RemissionLearnerMultiEdge_n(int dimension);
    
    inline void computeError();
    
    inline void linearizeOplusNumeric();
    inline void linearizeOplusJacobian();
    virtual void linearizeOplus();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
};

class RemissionCalibrationGridEdge : public g2o::BaseBinaryEdge<1, double, RemissionLearnerVertex, RemissionLearnerVertex>
{
  public:
    typedef g2o::BaseBinaryEdge<1, double, RemissionLearnerVertex, RemissionLearnerVertex> BaseClass;
    RemissionCalibrationGridEdge(RemissionLearnerVertex* v1, RemissionLearnerVertex* v2, double weight);
    inline void computeError();
    virtual void linearizeOplus();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class RemissionCalibrationEdge_2 : public g2o::BaseBinaryEdge<1, double, RemissionLearnerVertex, RemissionLearnerVertex>
{
  public:
    typedef g2o::BaseBinaryEdge<1, double, RemissionLearnerVertex, RemissionLearnerVertex> BaseClass;
    RemissionCalibrationEdge_2();
    inline void computeError();
    virtual void linearizeOplus();
    inline void linearizeOplusNumeric();
    inline void linearizeOplusJacobian();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class RemissionCalibrationEdgeMode2 : public g2o::BaseMultiEdge<1, double>
{
  public:
    typedef g2o::BaseMultiEdge<1, double> BaseClass;
    RemissionCalibrationEdgeMode2();
    inline void computeError();
    virtual void linearizeOplus();
    inline void linearizeOplusNumeric();
    inline void linearizeOplusJacobian();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class RemissionCalibrationEdgeMode3 : public g2o::BaseMultiEdge<1, double>
{
  public:
    typedef g2o::BaseMultiEdge<1, double> BaseClass;
    RemissionCalibrationEdgeMode3();
    inline void computeError();
    virtual void linearizeOplus();
    inline void linearizeOplusNumeric();
    inline void linearizeOplusJacobian();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

typedef std::map<int, std::vector<RemissionLearnerVertex*> > VerticesMode1;

struct CalibIndependentVertices {
  std::vector<RemissionLearnerVertex*> verticesRange,
                                       verticesIncidenceAngle,
                                       verticesIntensity;
};
typedef std::map<int, CalibIndependentVertices> VerticesMode2;

struct LaserGroupVertices {
  std::vector<RemissionLearnerVertex*> verticesRange;
  std::vector<RemissionLearnerVertex*> verticesIncidenceAngle;
  std::vector<RemissionLearnerVertex*> verticesIntensity;
  std::map<int, RemissionLearnerVertex*> verticesLaser;
};
typedef std::vector<LaserGroupVertices> VerticesMode3;

class G2oPostIterationAction : public g2o::HyperGraphAction {
  public:
    G2oPostIterationAction() : remissionCalibrationResult(NULL),
                               verticesMode1(NULL), verticesMode2(NULL), verticesMode3(NULL), verticesForMapCells(NULL),
                               mapCells(NULL) {}  //, qApplication(NULL), viewer(NULL) {}
    
    virtual g2o::HyperGraphAction* operator()(const g2o::HyperGraph* graph, g2o::HyperGraphAction::Parameters* parameters=NULL);
    
    RemissionCalibrationResult* remissionCalibrationResult;
    const VerticesMode1* verticesMode1;
    const VerticesMode2* verticesMode2;
    const VerticesMode3* verticesMode3;
    const std::vector<RemissionLearnerVertex*>* verticesForMapCells;
    RemissionCalibrationHelper::CollectedPointsMap* mapCells;
    //QApplication* qApplication;
    //Ais3dTools::ALUGLViewer* viewer;
};



#include "g2oGraphElements.hpp"

#endif
