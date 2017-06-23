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
	
	ListDigraph::NodeMap<bool> visitedNode(g);
	for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
		visitedNode[node] = false;
	}
	
	visitedNode[s] = true;
	
	double premioTotal = 0.0;
	
	ListDigraph::Node atual = s;
	
	ListDigraph::ArcMap<int> solInicial(g);
	for(ListDigraph::ArcIt arc(g); arc != INVALID; ++arc){
		solInicial[arc] = 0;
	}
	
	while(g.id(atual) != g.id(t)){
		double custo = INFINITY;
		ListDigraph::Arc arestaMenorCusto;
		bool achou = false;
		for (ListDigraph::OutArcIt aresta(g, atual); aresta != INVALID; ++aresta){
			if(g.id(g.target(aresta)) == g.id(t)){
				achou = true;
				arestaMenorCusto = aresta;
				custo = (-1)* INFINITY;
			}else{
				if(! visitedNode[g.target(aresta)]){
					if(prize[g.target(aresta)] - cost[aresta] < custo){
						achou = true;
						arestaMenorCusto = aresta;
						custo = cost[aresta];
					}
				}
			}
		}
		
		if(!achou){
			cout << "ERRO: Não achou caminho\n";
			return 0;
		}
		visitedNode[g.target(arestaMenorCusto)] = true;
		solInicial[arestaMenorCusto] = 1;
		atual = g.target(arestaMenorCusto);
		
		premioTotal += prize[g.target(arestaMenorCusto)] - cost[arestaMenorCusto];
	}
	
	// desconta premio de t
	premioTotal -= prize[t];
	
	// FIM DA HEURISTICA GULOSA
	
	clock_t posGuloso = clock();
	
	if( (double(posGuloso - begin) / CLOCKS_PER_SEC) > tMax ){
		UB = INFINITY;
		LB = premioTotal;
		return 2; // sol heuristica
	}
	
	// INICIO DO PLI
	
	try{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBLinExpr expr;
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		
		// o modelo ja comeca com o valor do premio da solucao gulosa como parametro de cutoff
		model.getEnv().set(GRB_DoubleParam_Cutoff, premioTotal - 1.0);
		
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
			
			//exprObjetivo += x[g.id(arc)]*(prize[g.target(arc)] - cost[arc]);
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
			UB = INFINITY;
			LB = premioTotal;
			return false;
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
		
		if (status == GRB_OPTIMAL){ // solucao otima
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
	return 0;
}
