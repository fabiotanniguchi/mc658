/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * EDITE ESTE ARQUIVO APENAS ONDE INDICADO
 * DIGITE SEU RA: 145980
 * SUBMETA SOMENTE ESTE ARQUIVO
 ******************************************************************************/

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tsp_bt_bnb.h"

bool bt(TSP_Data &tsp, int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BACKTRACKING.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	clock_t begin = clock();
	
	NodeIt nit(tsp.g);
	Node n = nit;
	
	// percorrer o grafo g fazendo backtracking a partir da sol gulosa
	greedy(tsp, maxTime);
	
	NodeBoolMap visitado(tsp.g);
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		visitado[o] = false;
	}
	
	clock_t posGuloso = clock();
	
	if( (double(posGuloso - begin) / CLOCKS_PER_SEC) > maxTime ){
		return false;
	}
	
	for (std::vector<Node>::iterator it = tsp.BestCircuit.begin() ; it != tsp.BestCircuit.end(); ++it){
		Node atual = *it;
		visitado[atual] = true;
		//cout << "Visitado: " << tsp.g.id(atual) << endl;
	}
	
	double minCost = tsp.BestCircuitValue;
	double costAux = minCost;
	
	for (std::vector<Node>::iterator it = tsp.BestCircuit.end() ; it != tsp.BestCircuit.begin() + 1; --it){
		Node atual = *it;
		Node antecessor = *(it - 1);
		Edge antToAtu = findEdge(tsp.g, antecessor, atual);
		costAux = costAux - tsp.weight[antToAtu];
		
		IncEdgeIt e(tsp.g, antecessor);
		for(; e != INVALID; ++e){
			Node adj = tsp.g.v(e);
			
			if(!visitado[adj]){
				int distance = std::distance( it, tsp.BestCircuit.end() );
				
				if(distance == 0){
					IncEdgeIt eAdj(tsp.g, adj);
					
					for(; eAdj != INVALID; ++eAdj){
						if(tsp.g.v(eAdj) == n){
							if(costAux + tsp.weight[eAdj] < minCost){
								minCost = costAux + tsp.weight[eAdj];
								tsp.BestCircuitValue = minCost;
								tsp.BestCircuit.pop_back();
								tsp.BestCircuit.push_back(adj);
							}
						}
					}
				}else{
					vector<Node> stack;
					vector<Edge> edgeStack;
					vector<Edge> solutionCandidateStack;
					
					IncEdgeIt eAdj(tsp.g, adj);
					
					for(; eAdj != INVALID; ++eAdj){
						Node destino = tsp.g.v(eAdj);
						if(!visitado[destino] && destino != n){
							stack.push_back(tsp.g.v(eAdj));
							
							Edge adjEdge = findEdge(tsp.g, adj, destino);
							edgeStack.push_back(adjEdge);
						}
					}
					
					while(stack.size() != 0){
						double internalCost = costAux;
						Node nodeFromStack = stack[stack.size() - 1];
						Edge edgeFromStack = edgeStack[edgeStack.size() - 1];
						stack.pop_back();
						edgeStack.pop_back();
						solutionCandidateStack.push_back(edgeFromStack);
						
						visitado[nodeFromStack] = true;
						
						IncEdgeIt eAdjNodeFromStack(tsp.g, nodeFromStack);
						bool achouVisitado = false;
						for(; eAdjNodeFromStack != INVALID; ++eAdjNodeFromStack){
							Node adjNode = tsp.g.v(eAdjNodeFromStack);
							
							if(!visitado[adjNode]){
								stack.push_back(adjNode);
								edgeStack.push_back(eAdjNodeFromStack);
								achouVisitado = true;
							}
						}
						
						if(!achouVisitado){ // chegou num caminho candidato a solucao do TSP
							for (std::vector<Edge>::iterator itEdges = solutionCandidateStack.begin() ; itEdges != solutionCandidateStack.end(); ++itEdges){
								internalCost += tsp.weight[*itEdges];
							}
							
							if(internalCost < minCost){
								Edge beginEdge = *(solutionCandidateStack.begin());
								Node intermediaryNode = tsp.g.u(beginEdge);
								while(tsp.BestCircuit[tsp.BestCircuit.size() - 1] != intermediaryNode && tsp.BestCircuit.size() > 0){
									tsp.BestCircuit.pop_back();
								}
								
								for (std::vector<Edge>::iterator itEdges = solutionCandidateStack.begin() ; itEdges != solutionCandidateStack.end() - 1; ++itEdges){
									Node nodeToAdd = tsp.g.v(*itEdges);
									tsp.BestCircuit.push_back(nodeToAdd);
								}
								
								tsp.BestCircuitValue = internalCost;
								minCost = internalCost;
							}
							
							visitado[nodeFromStack] = false;
							
							Edge last = edgeStack[edgeStack.size() - 1];
							Node uLast = tsp.g.u(last);
							
							while(tsp.g.v(solutionCandidateStack[solutionCandidateStack.size() - 1]) != uLast && solutionCandidateStack.size() > 0){
								// restaura situacao caso o backtracking tenha que voltar mais de uma aresta
								visitado[tsp.g.v(solutionCandidateStack[solutionCandidateStack.size() - 1])] = false;
								solutionCandidateStack.pop_back();
							}
						}
					}
				}
			}
			
			clock_t clockPosAnalise = clock();
	
			if( (double(clockPosAnalise - begin) / CLOCKS_PER_SEC) > maxTime ){
				return false;
			}
		}
		
		visitado[atual] = false;
	}
	
	return true;
}

bool bnb(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BRANCH AND BOUND.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	clock_t begin = clock();
	
	NodeIt nit(tsp.g);
	Node n = nit;
	
	// inicio da montagem da arvore de ramificacoes
	
	ListGraph ramtree; // arvore de ramificacoes
	ListGraph::NodeMap<Node> equivalencia(ramtree);
	
	Node raiz = ramtree.addNode();
	equivalencia[raiz] = n;
	
	vector<Node> stack;
	vector<Edge> edgeStack;
	vector<Edge> possibilityStack;
	
	NodeBoolMap visitado(tsp.g);
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		visitado[o] = false;
	}
	visitado[n] = true;
	
	IncEdgeIt eAdj(tsp.g, n);
	
	for(; eAdj != INVALID; ++eAdj){
		Node destino = tsp.g.v(eAdj);
		if(!visitado[destino] && destino != n){
			stack.push_back(tsp.g.v(eAdj));
			
			Edge adjEdge = findEdge(tsp.g, n, destino);
			edgeStack.push_back(adjEdge);
		}
	}

	while(stack.size() != 0){
		Node nodeFromStack = stack[stack.size() - 1];
		Edge edgeFromStack = edgeStack[edgeStack.size() - 1];
		stack.pop_back();
		edgeStack.pop_back();
		possibilityStack.push_back(edgeFromStack);
		
		visitado[nodeFromStack] = true;
		
		IncEdgeIt eAdjNodeFromStack(tsp.g, nodeFromStack);
		bool achouVisitado = false;
		for(; eAdjNodeFromStack != INVALID; ++eAdjNodeFromStack){
			Node adjNode = tsp.g.v(eAdjNodeFromStack);
			
			if(!visitado[adjNode]){
				stack.push_back(adjNode);
				edgeStack.push_back(eAdjNodeFromStack);
				achouVisitado = true;
			}
		}
		
		if(!achouVisitado){ // chegou num caminho candidato a solucao do TSP
			
			Node lastNode;
			Node novoNo = raiz;
			
			for (std::vector<Edge>::iterator itEdges = possibilityStack.begin() ; itEdges != possibilityStack.end(); ++itEdges){
				lastNode = novoNo;
				novoNo = ramtree.addNode();
				equivalencia[novoNo] = tsp.g.v(*itEdges);
				
				ramtree.addEdge(lastNode, novoNo);
			}
			
			visitado[nodeFromStack] = false;
			
			Edge last = edgeStack[edgeStack.size() - 1];
			Node uLast = tsp.g.u(last);
			
			while(tsp.g.v(possibilityStack[possibilityStack.size() - 1]) != uLast && possibilityStack.size() > 0){
				// restaura situacao caso o backtracking tenha que voltar mais de uma aresta
				visitado[tsp.g.v(possibilityStack[possibilityStack.size() - 1])] = false;
				possibilityStack.pop_back();
			}
		}
	}
	
	// termino da montagem da arvore de ramificacoes
	
	clock_t posMontagem = clock();
	
	if( (double(posMontagem - begin) / CLOCKS_PER_SEC) > maxTime ){
		return false;
	}
	
	// inicio da poda na arvore de ramificacoes partindo da solucao gulosa
	greedy(tsp, maxTime);
	
	IncEdgeIt e(ramtree, raiz);
	Edge minPathEdge = INVALID;
	double minCost = tsp.BestCircuitValue;
	
	for(; e != INVALID; ++e){
		Edge adjEdge = e;
		Node adj = ramtree.v(e);
		
		vector<Edge> candidateVector;
		double costAux = 0;
		Edge edgeAux = e;
		bool continuePath = true;
		bool achouMinCost = false;
		while(continuePath){
			candidateVector.push_back(edgeAux);
			Node u = equivalencia[ramtree.u(edgeAux)];
			Node v = equivalencia[ramtree.v(edgeAux)];
			Edge edgeEq = findEdge(tsp.g, u, v);
			costAux += tsp.weight[edgeEq];
			
			if(costAux >= minCost){
				continuePath =  false;
			}else{
				IncEdgeIt aux(ramtree, ramtree.v(edgeAux));
				if(aux == INVALID){ // chegou ao no folha
					continuePath = false;
					
					achouMinCost = true;
					
					tsp.BestCircuitValue = minCost;
					
					while(tsp.BestCircuit.size() > 0){
						tsp.BestCircuit.pop_back();
					}
					
					for (std::vector<Edge>::iterator itEdges = candidateVector.begin() ; itEdges != candidateVector.end(); ++itEdges){
						Node uEq = equivalencia[ramtree.u(*itEdges)];
						Node vEq = equivalencia[ramtree.v(*itEdges)];
						tsp.BestCircuit.push_back(uEq);
						
						if(itEdges == candidateVector.end() - 1){
							tsp.BestCircuit.push_back(vEq);
						}
					}
					
					if(minPathEdge != INVALID){
						ramtree.erase(minPathEdge);
						minPathEdge = adjEdge;
					}
				}else{
					edgeAux = aux;
				}
			}
		}
		
		if(!achouMinCost){
			ramtree.erase(adjEdge);
		}
		
		clock_t posLoopIt = clock();
	
		if( (double(posLoopIt - begin) / CLOCKS_PER_SEC) > maxTime ){
			return false;
		}
	}
	
	
	return true;
}

bool greedy(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * Algoritmo guloso para o TSP.
 ******************************************************************************/
{
	// Maps a boolean x to each edge of the graph g. Initialized as false.
	EdgeBoolMap x(tsp.g);
	for(EdgeIt e(tsp.g); e!=INVALID; ++e){  // We set every edge out of the the solution
		x[e] = false;
	}
	
	// NodeBoolMap y(tsp.g, false);  // Maps a boolean y to each vertex of the graph g. Initialized as false.
	NodeBoolMap y(tsp.g);
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		y[o] = false;
	}
	
	double cost = 0.0;
	int nedges = 0;
	
	//cout << endl;
	for(NodeIt o(tsp.g); o != INVALID; ++o){
		//cout << tsp.g.id(o) << ":" << y[o] << "  ";
	}
	//cout << endl;

	NodeIt nit(tsp.g);
	Node n = nit;
	Node f = nit;
	
	while(nedges != tsp.NNodes){
		// Put the vertex in the solution
		y[n] = true;
		tsp.BestCircuit.push_back(n);
		
		//cout << "n: " << tsp.g.id(n) << "  y[n]: " << y[n] << endl;
		
		double wmin = DBL_MAX;  // min weight
		IncEdgeIt emin = INVALID;  // min inc edge of n
		Node nmin = INVALID;
		
		IncEdgeIt e(tsp.g, n);
		Node op = INVALID;
		
		//cout << "wmin: " << wmin << endl;
		
		for(; e != INVALID; ++e){

			
			op = tsp.g.v(e);
			if( op == n ){
				op = tsp.g.u(e);
			}
			
			//cout << "   (" << tsp.g.id(tsp.g.u(e)) << ", " << tsp.g.id(tsp.g.v(e)) << ")  x: " << x[e] << "  c: " << tsp.weight[e] << " op: " << tsp.g.id(op) << " y[op] " << y[op] << endl;
			
			if( ! y[ op ] ){        // The expression in [] returns the "destin" vertex of edge e
				if( tsp.weight[e] < wmin ){
					wmin = tsp.weight[e];
					emin = e;
					nmin = op;
				}
			}
		}
		
		if( wmin < DBL_MAX ){  // If got some edge
			//cout << "wmin: " << wmin << endl;
			x[emin] = true;  // Puts the edge e in the solution: this data will be visible outside this function
			nedges++;
			cost += wmin;
			n = nmin;
			//cout << "new n: " << tsp.g.id(n) << endl;
		}
		else{
			//cout << "Error: could not found a minimum weight value." << endl;
			exit(1);
		}
		
		if( nedges == tsp.NNodes - 1 ){
			y[f] = false;
		}
		
		// cerr << "nedges: " << nedges << endl;
		// cerr << endl;
	}
	
	if( nedges > 0 ){
		tsp.BestCircuitValue = cost;
	}

	return false;
}
