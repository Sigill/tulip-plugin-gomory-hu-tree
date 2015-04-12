#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

#include <tulip/TulipPluginHeaders.h>
#include <tulip/ConnectedTest.h>
#include <tulip/tuliphash.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

#define CHECK_PROP_PROVIDED(PROP, STOR) \
	do { \
		if(!dataSet->get(PROP, STOR)) \
		throw std::runtime_error(std::string("No \"") + PROP + "\" property provided."); \
	} while(0)

namespace {
	const char * paramHelp[] = {
		HTML_HELP_OPEN() \
			HTML_HELP_DEF( "type", "IntegerProperty" ) \
			HTML_HELP_BODY() \
			"The capacity of the edges." \
			HTML_HELP_CLOSE()
	};
}

//#define DEBUG_PRINTF

using namespace std;
using namespace tlp;

bool compareNodes(const node a, const node b) {
	return a.id < b.id;
}

class IntegerMetaValueSumCalculator :public IntegerMinMaxProperty::MetaValueCalculator {
public:
	void computeMetaValue(AbstractProperty<IntegerType, IntegerType, NumericProperty>* metric, edge mE, Iterator<edge>* itE, Graph*) {
		int value = 0;

		while(itE->hasNext()) {
			edge e = itE->next();
			value += metric->getEdgeValue(e);
		}

		metric->setEdgeValue(mE, value);
	}
};

static IntegerMetaValueSumCalculator SumIntegerMetaCalculator;

typedef boost::adjacency_list_traits < boost::vecS, boost::vecS, boost::directedS > Traits;
typedef boost::adjacency_list<
		boost::vecS, boost::vecS, boost::directedS,

		boost::property < boost::vertex_index_t, long,
			boost::property < boost::vertex_index2_t, long,
			boost::property < boost::vertex_color_t, boost::default_color_type,
			boost::property < boost::vertex_distance_t, long,
			boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor > > > > >,

		boost::property < boost::edge_capacity_t, long,
			boost::property < boost::edge_residual_capacity_t, long,
			boost::property < boost::edge_reverse_t, Traits::edge_descriptor > > >
	> boost_directed_graph;

typedef Traits::vertex_descriptor vertex_t;
typedef Traits::edge_descriptor edge_t;

typedef boost::property_map <boost_directed_graph, boost::vertex_color_t >::type colormap_t;
typedef typename boost::property_traits<colormap_t>::value_type colorvalue_t;
typedef boost::color_traits<colorvalue_t> colortraits;

const string PLUGIN_NAME("GomoryHu");

class GomoryHu: public Algorithm {
private:
	tlp::IntegerProperty *m_Capacity, *m_Flow;

public:
	PLUGININFORMATIONS(PLUGIN_NAME, "Cyrille FAUCHEUX", "2015", "", "1.0", "")

		GomoryHu(PluginContext *context) : Algorithm(context)
	{
		addInParameter< IntegerProperty >("capacity", paramHelp[0], "capacity");
	}

	~GomoryHu() {}

	bool check(string &err) {
		try {
			if(dataSet == NULL)
				throw std::runtime_error("No dataset provided.");

			CHECK_PROP_PROVIDED("capacity", m_Capacity);

			m_Capacity->setMetaValueCalculator(&SumIntegerMetaCalculator);
		} catch (std::runtime_error &ex) {
			err.assign(ex.what());
			return false;
		}

		return true;
	}

	bool run() {
		Graph *G = graph->addCloneSubGraph("Original Graph"); // Backup

		// T = (Vt, Et)
		BooleanProperty *tmpSelection = graph->getLocalProperty<BooleanProperty>("tmpSelection");
		tmpSelection->setAllNodeValue(true);
		Graph *T = graph->addSubGraph(tmpSelection, "Gomory-Hu Tree");
		graph->delLocalProperty("tmpSelection");

		// Set Vt = {Vg} and Et = {}
		{
			std::set<node> nodes;
			node n;
			forEach(n, T->getNodes())
				nodes.insert(n);

			T->createMetaNode(nodes);
			// TODO May be shorter to use the other version of createMetaNode
		}

#ifdef DEBUG_PRINTF
		int i = 0;
#endif
		while(true) {
#ifdef DEBUG_PRINTF
			std::cout << i++ << std::endl;
#endif
			// Choose some X in Vt with |X| > 1
			node X;
			Iterator<node> *itNodes = T->getNodes();
			while(itNodes->hasNext()) {
				node n = itNodes->next();
				if(T->getNodeMetaInfo(n)->numberOfNodes() > 1) {
					X = n;
					break;
				}
			} delete itNodes;

			// No more cut to compute
			if(!X.isValid())
				break;

#ifdef DEBUG_PRINTF
			node n;
			std::cout << "X = {";
			forEach(n, T->getNodeMetaInfo(X)->getNodes()) {
				std::cout << n.id << " ";
			}
			std::cout << "}" << std::endl;
#endif

			set<node> SC_keys;
			Graph *Gprime;
				Gprime = buildGprime(G, T, X, &SC_keys);

			itNodes = Gprime->getNodes();
			node n1 = itNodes->next();
			node n2 = itNodes->next();
			delete itNodes;

			boost_directed_graph g(Gprime->numberOfNodes());

			boost::property_map <boost_directed_graph, boost::vertex_index2_t >::type                   tulip_id = boost::get(boost::vertex_index2, g);
			boost::property_map <boost_directed_graph, boost::vertex_color_t >::type                       color = boost::get(boost::vertex_color, g);

			boost::property_map <boost_directed_graph, boost::edge_capacity_t >::type                   capacity = boost::get(boost::edge_capacity, g);
			boost::property_map <boost_directed_graph, boost::edge_residual_capacity_t >::type residual_capacity = boost::get(boost::edge_residual_capacity, g);
			boost::property_map <boost_directed_graph, boost::edge_reverse_t >::type                         rev = boost::get(boost::edge_reverse, g);

			vertex_t s, t;
			edge_t u, v;

			TLP_HASH_MAP<unsigned int, vertex_t> tulip2boost_map;
			node n;
			forEach(n, Gprime->getNodes()) {
				s = boost::add_vertex(g);
				boost::put(tulip_id, s, (long)(n.id));
				tulip2boost_map.insert(std::make_pair<unsigned int, vertex_t>(n.id, s));
			}

			bool inserted;
			int c;
			edge e;
			forEach(e, Gprime->getEdges()) {
				s = tulip2boost_map.find(Gprime->source(e))->second;
				t = tulip2boost_map.find(Gprime->target(e))->second;

				boost::tie(u, inserted) = boost::add_edge(s, t, g);
				boost::tie(v, inserted) = boost::add_edge(t, s, g);

				c = m_Capacity->getEdgeValue(e);
				boost::put(capacity, u, c);
				boost::put(capacity, v, c);

				rev[u] = v;
				rev[v] = u;
			}

			std::vector<long> distance(num_vertices(g));
			long flow = boost::boykov_kolmogorov_max_flow(g, tulip2boost_map.find(n1.id)->second, tulip2boost_map.find(n2.id)->second);

			MutableContainer<bool> markedNodes;
			boost::graph_traits < boost_directed_graph >::vertex_iterator u_iter, u_end;
			for (tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter) {
				if(get(color, *u_iter) == colortraits::black())
					markedNodes.set(node(boost::get(tulip_id, *u_iter)), true);
			}

#ifdef DEBUG_PRINTF
			std::cout << "Cut " << n1.id << "/" << n2.id << " = " << flow << std::endl;

			forEach(n, Gprime->getNodes()) {
				if(markedNodes.get(n))
					std::cout << n.id << " ";
			}
			std::cout << std::endl;

			forEach(n, Gprime->getNodes()) {
				if(!markedNodes.get(n))
					std::cout << n.id << " ";
			}
			std::cout << std::endl;
#endif

			updateT(&markedNodes, flow, G, Gprime, &SC_keys, T, X);
		}

		StableIterator<node> itMeta1(T->getNodes(), T->numberOfNodes(), true),
		                     itMeta2(T->getNodes(), T->numberOfNodes(), true);
		while(itMeta1.hasNext()) {
			T->addNode(graph->getNodeMetaInfo(itMeta1.next())->getOneNode());
		}

		Iterator<edge> *ite = T->getEdges();
		while(ite->hasNext()) {
			edge e = ite->next();
			node s = T->source(e),
			     t = T->target(e);

			T->setEnds(e, graph->getNodeMetaInfo(s)->getOneNode(), graph->getNodeMetaInfo(t)->getOneNode());
		} delete ite;

		while(itMeta2.hasNext()) {
			node n = itMeta2.next();
			T->delNode(n);
			graph->delAllSubGraphs(graph->getNodeMetaInfo(n));
		}

		return true;
	}

	Graph* buildGprime(Graph *G, Graph *T, node X, set<node> *SC_keys) {
		Graph *TminusX = T->addCloneSubGraph("T\\X");
		TminusX->delNode(X);

		std::vector< std::set<node> > CC;
		ConnectedTest::computeConnectedComponents(TminusX, CC);

		T->delAllSubGraphs(TminusX);

		Graph *Gprime = G->addCloneSubGraph("G'");

		// Contracting the connected components from T\X in G'
		for(std::vector< std::set<node> >::iterator cc = CC.begin() ; cc != CC.end(); ++cc) {
			set<node> nodes;
			for(std::set<node>::iterator mn = cc->begin() ; mn != cc->end(); ++mn) {
				node n;
				forEach(n, graph->getNodeMetaInfo(*mn)->getNodes())
					nodes.insert(n);
			}

			node n = Gprime->createMetaNode(nodes, true, false);
			SC_keys->insert(n);
		}

		// Merging the resulting edges in G'
		node n;
		forEach(n, graph->getNodeMetaInfo(X)->getNodes()) {
			for(set<node>::iterator mn = SC_keys->begin() ; mn != SC_keys->end(); ++mn) {
				double cost = 0;
				std::vector<edge> edges = Gprime->getEdges(n, *mn, false);
				for(std::vector< edge >::iterator it = edges.begin() ; it != edges.end(); ++it) {
					cost += m_Capacity->getEdgeValue(*it);
					Gprime->delEdge(*it);
				}

				if(cost > 0) {
					edge e = Gprime->addEdge(n, *mn);
					m_Capacity->setEdgeValue(e, cost);
				}
			}
		}

		return Gprime;
	}

	void updateT(MutableContainer<bool> *markedNodes, double cutCost, Graph *G, Graph* Gprime, set<node> *SC_keys, Graph *T, node X) {
		set<node> Xnodes;
		node n;
		forEach(n, graph->getNodeMetaInfo(X)->getNodes())
			Xnodes.insert(n);

		set<node> Aprime, Bprime;
		forEach(n, Gprime->getNodes()) {
			if(markedNodes->get(n)) {
				Aprime.insert(n);
			} else {
				Bprime.insert(n);
			}
		}

		set<node> A;
		for(set<node>::iterator n = Xnodes.begin() ; n != Xnodes.end(); ++n) { // A' inter X
			if(Aprime.find(*n) != Aprime.end()) {
				A.insert(*n);
			}
		}
		for(set<node>::iterator n = Aprime.begin() ; n != Aprime.end(); ++n) { // + Union SC in {A' inter S}
			if(SC_keys->find(*n) != SC_keys->end()) {
				node innerNode;
				forEach(innerNode, graph->getNodeMetaInfo(*n)->getNodes()) {
					A.insert(innerNode);
				}
			}
		}

		set<node> B;
		for(set<node>::iterator n = Xnodes.begin() ; n != Xnodes.end(); ++n) {
			if(Bprime.find(*n) != Bprime.end()) {
				B.insert(*n);
			}
		}
		for(set<node>::iterator n = Bprime.begin() ; n != Bprime.end(); ++n) {
			if(SC_keys->find(*n) != SC_keys->end()) {
				node innerNode;
				forEach(innerNode, graph->getNodeMetaInfo(*n)->getNodes()) {
					B.insert(innerNode);
				}
			}
		}

		for(set<node>::iterator n = SC_keys->begin() ; n != SC_keys->end(); ++n) {
			Graph *s = G->getNodeMetaInfo(*n);
			graph->delNode(*n, true);
			G->delAllSubGraphs(s);
		}

		G->delAllSubGraphs(Gprime);

		set<node> AinterX, BinterX;
		for(set<node>::iterator n = Xnodes.begin() ; n != Xnodes.end(); ++n) {
			if(Aprime.find(*n) != Aprime.end()) {
				AinterX.insert(*n);
			}
			if(Bprime.find(*n) != Bprime.end()) {
				BinterX.insert(*n);
			}
		}

		for(set<node>::iterator n = AinterX.begin() ; n != AinterX.end(); ++n) {
			T->addNode(*n);
		}

		node AinterX_node = T->createMetaNode(AinterX, false, false);
		graph->delEdges(graph->getInOutEdges(AinterX_node), true); // Those edges should not have been added

		for(set<node>::iterator n = BinterX.begin() ; n != BinterX.end(); ++n) {
			T->addNode(*n);
		}

		node BinterX_node = T->createMetaNode(BinterX, false, false);
		graph->delEdges(graph->getInOutEdges(BinterX_node), true); // Those edges should not have been added

		Iterator<edge> *pit = T->getInOutEdges(X);
		while(pit->hasNext()) {
			edge e = pit->next();
			node Y = graph->opposite(e, X);
			set<node> YNodes;
			if(T->isMetaNode(Y)) {
				node n;
				forEach(n, graph->getNodeMetaInfo(Y)->getNodes())
					YNodes.insert(n);
			} else {
				YNodes.insert(Y);
			}

			bool found = true;
			for(set<node>::iterator n = YNodes.begin() ; n != YNodes.end(); ++n) {
				if(A.find(*n) == A.end()) {
					found = false;
					break;
				}
			}

			node newX = found ? AinterX_node : BinterX_node;
			edge eprime = T->addEdge(newX, Y);
			m_Capacity->setEdgeValue(eprime, m_Capacity->getEdgeValue(e));

		} delete pit;

		m_Capacity->setEdgeValue(T->addEdge(AinterX_node, BinterX_node), cutCost);

		Graph *Xsubgraph = graph->getNodeMetaInfo(X);
		graph->delNode(X, true);
		graph->delAllSubGraphs(Xsubgraph);
	}

};

PLUGIN(GomoryHu);
