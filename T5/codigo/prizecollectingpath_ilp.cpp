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
int prize_collecting_st_path_pli(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double>& cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &UB, double &LB, int tMax){
	
	clock_t begin = clock();
	
	// pega valor de cutoff da heuristica
	double cutoffValue = prize_collecting_st_path_heuristic(g, prize, cost, s, t, path, UB, LB, tMax);
	
	clock_t posGuloso = clock();
	
	if( (double(posGuloso - begin) / CLOCKS_PER_SEC) > tMax ){
		UB = DBL_MAX;
		LB = cutoffValue;
		return 2; // sol heuristica
	}
	
	path.clear();
	
	// INICIO DO PLI
	
	try{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBLinExpr expr;
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		
		// o modelo ja comeca com o valor do premio da solucao gulosa como parametro de cutoff
		model.getEnv().set(GRB_DoubleParam_Cutoff, cutoffValue - 1.0);
		
		// criando variaveis
		vector<GRBVar> x(countArcs(g)); // variavel que indica se a aresta eh usada
		for(ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
			x[g.id(arc)] = model.addVar(0, 1, 0, GRB_BINARY);
		}
		
		model.update();
		
		// definindo objetivo
		GRBLinExpr exprObjetivo;
		for(ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
			exprObjetivo += x[g.id(arc)]*prize[g.target(arc)] - x[g.id(arc)]*cost[arc];
		}
		exprObjetivo = exprObjetivo - prize[t];
		
		model.setObjective(exprObjetivo, GRB_MAXIMIZE);
		
		model.update();
		
		// definindo restricoes
		for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
			GRBLinExpr restr1;
			GRBLinExpr restr2;
			GRBLinExpr restr3;
			for(ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc){
				ListDigraph::Arc original = arc;
				restr1 += x[g.id(original)];
				restr2 += x[g.id(original)];
			}
			for(ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc){
				ListDigraph::Arc original = arc;
				restr1 -= x[g.id(original)];
				restr3 += x[g.id(original)];
			}
			
			if(g.id(node) != g.id(s) && g.id(node) != g.id(t)){
				// grau de entrada deve ser igual ao grau de saida
				model.addConstr(restr1 == 0.0);
				
				// grau de entrada deve ser no maximo 1
				model.addConstr(restr2 <= 1.0);
				
				// grau de saida deve ser no maximo 1
				model.addConstr(restr3 <= 1.0);
			}else{
				if(g.id(node) == g.id(s)){ // no de partida
					// grau de entrada deve ser zero
					model.addConstr(restr2 == 0.0);
					
					// grau de saida deve ser igual a 1
					model.addConstr(restr3 == 1.0);
					
				}else{ // no de chegada do caminho
					// grau de entrada deve ser igual a 1
					model.addConstr(restr2 == 1.0);
					
					// grau de saida deve ser zero
					model.addConstr(restr3 == 0.0);
				}
			}
		}
		
		clock_t posFormulacao = clock();
		if( (double(posFormulacao - begin) / CLOCKS_PER_SEC) > tMax ){
			UB = DBL_MAX;
			LB = cutoffValue;
			return 2; // sol heuristica
		}else{
			double tUsed = (posFormulacao - begin) / CLOCKS_PER_SEC;
			// restrinjo o Gurobi com o tempo restante de execucao deste algoritmo
			model.getEnv().set(GRB_DoubleParam_TimeLimit, tMax - tUsed);
		}
		
		model.update();
		model.optimize(); // roda o solver
		int status = model.get(GRB_IntAttr_Status);
		
		ListDigraph::Node atual = s;
		path.push_back(s);
		
		double total = 0.0;
		
		// percorre de s a t passando pelos arcos escolhidos pelo solver
		// veja que isto remove os possiveis ciclos isolados gerados pelo PLI
		while(g.id(atual) != g.id(t)){
			for(ListDigraph::OutArcIt arc(g, atual); arc != INVALID; ++arc){
				ListDigraph::Arc theArc = arc;

				double originalValue = x[g.id(theArc)].get(GRB_DoubleAttr_X);
				long value = lround(originalValue);
				if(value == 1){
					path.push_back(g.target(theArc));
					atual = g.target(theArc);
					total += prize[g.target(theArc)] - cost[theArc];
				}
			}
		}
		total -= prize[t];
		
		clock_t end = clock();
		double tempo = (double)(end - begin) / CLOCKS_PER_SEC;
		cout.precision(4);
		cout << "Tempo: " << fixed << tempo << endl;
		cout << "Premio obtido: " << fixed << total << endl;
		
		if (status == GRB_OPTIMAL){
			LB = total;
			UB = total;
			return 1; // sol otima
		}else{
			LB = total;
			UB = model.get(GRB_DoubleAttr_MaxBound);;
			return 2; // sol heuristica
		}
		
	}catch(GRBException e) {
		cerr << "Nao foi possivel resolver o PLI." << endl;
		cerr << "Codigo de erro = " << e.getErrorCode() << endl;
		cerr << e.getMessage();
	}
	
	return 0;
}


///
//
// Heuristic function
///
int prize_collecting_st_path_heuristic(ListDigraph& g, ListDigraph::NodeMap<double>& prize, ListDigraph::ArcMap<double> &cost, ListDigraph::Node s, ListDigraph::Node t, std::vector<ListDigraph::Node> &path, double &UB, double &LB, int tMax){
	
	// HEURISTICA GULOSA
	
	ListDigraph::NodeMap<bool> visitedNode(g);
	for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
		visitedNode[node] = false;
	}
	
	visitedNode[s] = true;
	
	double premioTotal = 0.0;
	
	ListDigraph::Node atual = s;
	path.push_back(s);
	
	while(g.id(atual) != g.id(t)){
		double premio = DBL_MIN;
		ListDigraph::Arc arestaMaiorPremio;
		bool achou = false;
		for (ListDigraph::OutArcIt aresta(g, atual); aresta != INVALID; ++aresta){
			if(g.id(g.target(aresta)) == g.id(t)){
				achou = true;
				arestaMaiorPremio = aresta;
				premio = DBL_MAX;
			}else{
				if(! visitedNode[g.target(aresta)]){
					if(prize[g.target(aresta)] - cost[aresta] > premio){
						achou = true;
						arestaMaiorPremio = aresta;
						premio = prize[g.target(aresta)] - cost[aresta];
					}
				}
			}
		}
		
		if(!achou){
			cerr << "Nao ha caminho\n";
			return 0;
		}else{
			visitedNode[g.target(arestaMaiorPremio)] = true;
			path.push_back(g.target(arestaMaiorPremio));
			atual = g.target(arestaMaiorPremio);
			premioTotal += prize[g.target(arestaMaiorPremio)] - cost[arestaMaiorPremio];
		}
	}
	
	// desconta premio de t
	premioTotal -= prize[t];
	
	return premioTotal;
}
