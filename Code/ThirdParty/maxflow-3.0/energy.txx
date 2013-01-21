/*
 * energy.txx
 *
 *  Created on: Jan 21, 2013
 *      Author: oesteban
 */
#include "energy.h"

#ifndef ENERGY_TXX_
#define ENERGY_TXX_

template <typename captype, typename tcaptype, typename flowtype>
inline Energy<captype,tcaptype,flowtype>::Energy(int var_num_max, int edge_num_max, void (*err_function)(const char *)) : Graph<captype,tcaptype,flowtype>(var_num_max, edge_num_max, err_function)
{
	Econst = 0;
	error_function = err_function;
}

template <typename captype, typename tcaptype, typename flowtype>
inline Energy<captype,tcaptype,flowtype>::~Energy() {}

template <typename captype, typename tcaptype, typename flowtype>
inline typename Energy<captype,tcaptype,flowtype>::Var Energy<captype,tcaptype,flowtype>::add_variable(int num)
{	return GraphT::add_node(num); }

template <typename captype, typename tcaptype, typename flowtype>
inline void Energy<captype,tcaptype,flowtype>::add_constant(Value A) { Econst += A; }

template <typename captype, typename tcaptype, typename flowtype>
inline void Energy<captype,tcaptype,flowtype>::add_term1(Var x,
                              Value A, Value B)
{
	GraphT::add_tweights(x, B, A);
}

template <typename captype, typename tcaptype, typename flowtype>
inline void Energy<captype,tcaptype,flowtype>::add_term2(Var x, Var y,
                              Value A, Value B,
                              Value C, Value D)
{
	/*
	   E = A A  +  0   B-A
	       D D     C-D 0
	   Add edges for the first term
	*/
	GraphT::add_tweights(x, D, A);
	B -= A; C -= D;

	/* now need to represent
	   0 B
	   C 0
	*/

	assert(B + C >= 0); /* check regularity */
	if (B < 0)
	{
		/* Write it as
		   B B  +  -B 0  +  0   0
		   0 0     -B 0     B+C 0
		*/
		GraphT::add_tweights(x, 0, B); /* first term */
		GraphT::add_tweights(y, 0, -B); /* second term */
		GraphT::add_edge(x, y, 0, B+C); /* third term */
	}
	else if (C < 0)
	{
		/* Write it as
		   -C -C  +  C 0  +  0 B+C
		    0  0     C 0     0 0
		*/
		GraphT::add_tweights(x, 0, -C); /* first term */
		GraphT::add_tweights(y, 0, C); /* second term */
		GraphT::add_edge(x, y, B+C, 0); /* third term */
	}
	else /* B >= 0, C >= 0 */
	{
		GraphT::add_edge(x, y, B, C);
	}
}

template <typename captype, typename tcaptype, typename flowtype>
inline void Energy<captype,tcaptype,flowtype>::add_term3(Var x, Var y, Var z,
                              Value E000, Value E001,
                              Value E010, Value E011,
                              Value E100, Value E101,
                              Value E110, Value E111)
{
	register Value pi = (E000 + E011 + E101 + E110) - (E100 + E010 + E001 + E111);
	register Value delta;
	register Var u;

	if (pi >= 0)
	{
		Econst += E111 - (E011 + E101 + E110);

		add_tweights(x, E101, E001);
		add_tweights(y, E110, E100);
		add_tweights(z, E011, E010);

		delta = (E010 + E001) - (E000 + E011); /* -pi(E[x=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(y, z, delta, 0);

		delta = (E100 + E001) - (E000 + E101); /* -pi(E[y=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(z, x, delta, 0);

		delta = (E100 + E010) - (E000 + E110); /* -pi(E[z=0]) */
		assert(delta >= 0); /* check regularity */
		add_edge(x, y, delta, 0);

		if (pi > 0)
		{
			u = add_variable();
			add_edge(x, u, pi, 0);
			add_edge(y, u, pi, 0);
			add_edge(z, u, pi, 0);
			add_tweights(u, 0, pi);
		}
	}
	else
	{
		Econst += E000 - (E100 + E010 + E001);

		add_tweights(x, E110, E010);
		add_tweights(y, E011, E001);
		add_tweights(z, E101, E100);

		delta = (E110 + E101) - (E100 + E111); /* -pi(E[x=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(z, y, delta, 0);

		delta = (E110 + E011) - (E010 + E111); /* -pi(E[y=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(x, z, delta, 0);

		delta = (E101 + E011) - (E001 + E111); /* -pi(E[z=1]) */
		assert(delta >= 0); /* check regularity */
		add_edge(y, x, delta, 0);

		u = add_variable();
		add_edge(u, x, -pi, 0);
		add_edge(u, y, -pi, 0);
		add_edge(u, z, -pi, 0);
		add_tweights(u, -pi, 0);
	}
}

template <typename captype, typename tcaptype, typename flowtype>
inline typename Energy<captype,tcaptype,flowtype>::TotalValue Energy<captype,tcaptype,flowtype>::minimize() {
return Econst + GraphT::maxflow(); }

template <typename captype, typename tcaptype, typename flowtype>
inline int Energy<captype,tcaptype,flowtype>::get_var(Var x) { return (int) GraphT::what_segment(x); }



#endif /* ENERGY_TXX_ */
