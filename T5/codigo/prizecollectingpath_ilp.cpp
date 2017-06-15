/*******************************************************
 * MC658 - Projeto e Analise de Algoritmo III - 1s2017
 * Prof: Flavio Keidi Miyazawa
 * PED: Edson Ticona Zegarra
 ******************************************************/
#include "prizecollectingpath.h"

///Preencher aqui para facilitar a correcao. 
// Nome1: Fabio Takahashi Tanniguchi
// RA1: 145980
// Nome2: -------------------------
// RA2: ---------------------------

///
// PLI function
///
int prize_collecting_st_path_pli(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBLinExpr expr;
	model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
	env.set(GRB_DoubleParam_TimeLimit, tMax);
	
	ListDigraph::NodeMap<bool> visitedNode(g);
	for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
		visitedNode[node] = false;
	}
	
	visitedNode[s] = true;
	
	vector<ListDigraph::Node> solInicial;
	solInicial.push_back(s);
	
	ListDigraph::Node atual = s;
	
	while(g.id(atual) != g.id(t)){
		double custo = INFINITY;
		ListDigraph::Arc arestaMenorCusto;
		bool achou = false;
		for (ListDigraph::OutArcIt aresta(g, atual); aresta != INVALID; ++aresta){
			if(! visitedNode[g.target(aresta)]){
				if(cost[aresta]-prize[g.target(aresta)] < custo){
					achou = true;
					arestaMenorCusto = aresta;
				}
			}
		}
		
		if(!achou){
			cout << "Não achou caminho\n";
			return 0;
		}
		visitedNode[g.target(arestaMenorCusto)] = true;
		solInicial.push_back(g.target(arestaMenorCusto));
	}
	
	return 0;
}


///
//
// Heuristic function
///
int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &LB, double &UB, int tMax){
	return 0;
}
