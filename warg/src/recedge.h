#ifndef RECEDGE_H
#define RECEDGE_H

namespace weakarg
{

/**
    @brief Recombinant edge
*/

class RecEdge
{
    friend class RecTree;
protected:
    double timeFrom,timeTo;///< times of the edge in coalescent time
    unsigned int gstart,gend;///< positions effected by the edge on the genome
    unsigned int edgeFrom,edgeTo;///< Edges on the coalescent tree affected
public:
    RecEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,unsigned int edgefrom,unsigned int edgeto);///< Creates a specified edge
    void MoveEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,unsigned int edgefrom,unsigned int edgeto);///< Moves the edge as specified
    ~RecEdge();
	RecEdge(const RecEdge& r) :
		timeFrom(r.timeFrom),
		timeTo(r.timeTo),
		gstart(r.gstart),
		gend(r.gend),
		edgeFrom(r.edgeFrom),
		edgeTo(r.edgeTo)
	{}

	 inline double getTimeFrom() const
    {
        return timeFrom;
    }
    ;///<Returns the time an edge starts at
    inline double getTimeTo() const
    {
        return timeTo;
    }
    ;///<Returns the time an edge arrives at
    inline unsigned int getStart() const
    {
        return gstart;
    }
    ;///<Returns the genome start point
    inline unsigned int getEnd() const
    {
        return gend;
    }
    ;///<Returns the genome end point
    inline int getEdgeFrom() const
    {
        return edgeFrom;
    }
    ;///<Returns the edge that material in this edge is from
    inline int getEdgeTo() const
    {
        return edgeTo;
    }
    ;///<Returns the edge that material is this edge is placed in
    inline void setTimeFrom(double t)
    {
        timeFrom=t;
    }
    ;///<Sets the time an edge leaves from
    inline void setTimeTo(double t)
    {
        timeTo=t;
    }
    ;///<Sets the time an edge arrives at
    inline void setEdgeFrom(unsigned int e)
    {
        edgeFrom=e;
    }
    ;///<Sets the edge an edge leaves from
    inline void setEdgeTo(unsigned int e)
    {
        edgeTo=e;
    }
    ;///<Sets the edge an edge arrives at
    inline void setTimeFrom(double t,unsigned int e)
    {
        timeFrom=t;
        edgeFrom=e;
    }
    ;///<Sets the time and edge an edge leaves from
    inline void setTimeTo(double t,unsigned int e)
    {
        timeTo=t;
        edgeTo=e;
    }
    ;///<Sets the time and edge an edge arrives at
    inline bool affectsSite(unsigned int i) const
    {
        if((i<gend)&&(i>=gstart))
            return(1);
        return(0);
    }///< Returns whether the edge affects a given site
};

} // end namespace weakarg
#endif
