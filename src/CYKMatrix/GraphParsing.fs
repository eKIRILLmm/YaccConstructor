﻿module GraphParsing
    open Util
    open TPMatrices
    open System
    open System.Collections.Generic
    open Yard.Core
    open Yard.Core.IL
    open Yard.Core.IL.Production
    open Yard.Core.IL.Definition
    open Yard.Core.Helpers
    open Conversions.TransformAux

    module Graph = 
        type T(_graph: Dictionary<int * int, ResizeArray<int>>, _numberOfVertices: int) =
            member this.graph = _graph

            member this.numberOfVertices = _numberOfVertices

            member this.IsCorrectVertex (v: int) =
                (1 <= v) && (v <= this.numberOfVertices)

            member this.AddEdge (v1: int, v2: int, label: int) =
                assert (this.IsCorrectVertex(v1) && this.IsCorrectVertex(v2))
                if not <| this.graph.ContainsKey (v1, v2)
                then
                    this.graph.Add((v1, v2), new ResizeArray<int>())
                if not <| this.graph.[(v1,v2)].Contains label
                then
                    this.graph.[(v1,v2)].Add(label)

            member this.GetEdges () =
                this.graph.Keys

            member this.GetLabels (v1: int, v2: int) =
                assert (this.IsCorrectVertex(v1) && this.IsCorrectVertex(v2))
                this.graph.[(v1, v2)]


    type ParsingMatrix = Dictionary<NonTerminal, ProbabilityMatrix.T>

    let initParsingMatrix (graph:Graph.T)
                  (allRules: RulesHolder)
                  nonterminals =
        let unionWithOneCell (arr: float []) (v1: int) (v2: int) (matrixSize: int) =
            Array.init (matrixSize * matrixSize) (fun i -> if arr.[i] > 0.0
                                                           then 1.0
                                                           else 
                                                               let row = i / matrixSize
                                                               let col = i - row * matrixSize
                                                               if (v1 - 1 = row) && (v2 - 1 = col)
                                                               then 1.0
                                                               else 0.0)
        let parsingMatrix = new ParsingMatrix ()
        do 
            (
                nonterminals 
                |> Seq.map (fun x -> x, ProbabilityMatrix.empty (graph.numberOfVertices))
            )
            |> Seq.iter parsingMatrix.Add

        for (v1, v2) in graph.GetEdges() do
            let labels = graph.GetLabels(v1, v2)
            for label in labels do
                if allRules.IsSimpleTail label
                then
                    let simpleNonterminals = allRules.HeadsBySimpleTail label
                    for (simpleNonterminal, _) in simpleNonterminals do
                        let arr = parsingMatrix.[simpleNonterminal].GetSubArray id false parsingMatrix.[simpleNonterminal].WholeMatrix
                        let updatedArr = unionWithOneCell arr v1 v2 graph.numberOfVertices
                        let updatedMatrix = ProbabilityMatrix.create parsingMatrix.[simpleNonterminal].Size updatedArr
                        parsingMatrix.Remove(simpleNonterminal) |> ignore
                        parsingMatrix.Add(simpleNonterminal, updatedMatrix)
        
        parsingMatrix

    let naiveSquareMatrix (matrix: ParsingMatrix) (allRules: RulesHolder) isChanged =
        let unionArrays (arr1: float []) (arr2: float []) (size: int) =
            let modificator i =
                if arr1.[i] <> arr2.[i]
                then isChanged := true

                if arr1.[i] > 0.0
                then 1.0
                else arr2.[i]
            Array.init size modificator

        let multArrays (from1: Probability.InnerType.T []) (from2: Probability.InnerType.T []) matricesSize actualColCount =        
                let calculateCell (n, i, j) = 
                    let skipMatrices = n * matricesSize * matricesSize
                    let skipRows = skipMatrices + i * matricesSize
                    let skipColumns = skipMatrices + j * matricesSize                
                    Array.fold2 (fun b v1 v2 -> Probability.innerSumm b <| Probability.innerMult v1 v2)
                                Probability.InnerType.zero
                                from1.[skipRows..skipRows + matricesSize - 1] 
                                from2.[skipColumns..skipColumns + matricesSize - 1]

                // todo: recurr
                let getNIJ x = 
                    let n = x / (matricesSize * matricesSize)
                    let i = (x - n * matricesSize * matricesSize) / matricesSize
                    let j = x - n * matricesSize * matricesSize - i * matricesSize
                    n, i, j
                Array.init (matricesSize * actualColCount) (fun x -> calculateCell <| getNIJ x)

        let nontermPairs = allRules.ComplexTails
        for (nt1, nt2) in nontermPairs do
            let size = matrix.[nt1].Size
            let arr1 = matrix.[nt1].GetSubArray id false matrix.[nt1].WholeMatrix
            let arr2 = matrix.[nt2].GetSubArray id true matrix.[nt2].WholeMatrix
            let resultArray = multArrays arr1 arr2 size size
            let resultMatix = ProbabilityMatrix.create size resultArray

            for (nonTerm, _) in allRules.HeadsByComplexTail (nt1, nt2) do
                let curArray = matrix.[nonTerm].GetSubArray id false matrix.[nonTerm].WholeMatrix
                let updatedArray = unionArrays (resultMatix.GetSubArray id false resultMatix.WholeMatrix) curArray curArray.Length
                let updatedMatrix = ProbabilityMatrix.create matrix.[nonTerm].Size updatedArray
                matrix.Remove(nonTerm) |> ignore
                matrix.Add(nonTerm, updatedMatrix)


    let recognizeGraph (graph:Graph.T)
                  squareMatrix
                  (allRules: RulesHolder)
                  nonterminals
                  S =
        let parsingMatrix = initParsingMatrix graph allRules nonterminals
        let isChanged = ref true
        let mutable multCount = 0

        while !isChanged do
            isChanged := false
            squareMatrix parsingMatrix allRules isChanged
            multCount <- multCount + 1

        (parsingMatrix.[S], multCount)


    let graphParse (graph:Graph.T)
                  squareMatrix
                  (loadIL:t<Source.t, Source.t>) =
        let grammar = loadIL.grammar
        
        let S = ref (NonTerminal "")
        let nonterminals = new ResizeArray<NonTerminal>()
        let crl = new Dictionary<NonTerminal * NonTerminal, ResizeArray<NonTerminal*Probability.T>>()
        let srl = new Dictionary<int, ResizeArray<NonTerminal*Probability.T>>()
        let crl_result = new Dictionary<NonTerminal * NonTerminal, (NonTerminal * Probability.T) list>()
        let srl_result = new Dictionary<int, (NonTerminal * Probability.T) list>()
        let erl_result: NonTerminal list = []

        let probOne = Probability.create 1.0

        for module' in grammar do
            for r in module'.rules do
                let nonterm = NonTerminal <| Source.toString r.name
                if not <| nonterminals.Contains nonterm
                then
                    nonterminals.Add nonterm
                    if r.isStart
                    then
                        S := nonterm

                match r.body with
                | PSeq([elem],_,_) ->
                    match elem.rule with
                    | PToken src ->
                        let token = (Source.toString src).[0] |> int
                        if not <| srl.ContainsKey token
                        then
                            srl.Add(token, new ResizeArray<NonTerminal*Probability.T>())
                        if not <| srl.[token].Contains (nonterm, probOne)
                        then
                            srl.[token].Add (nonterm, probOne)
                    | _ ->
                        assert false
                        
                | PSeq([e1; e2],_,_) ->
                    match e1.rule, e2.rule with 
                    | PRef (name1, _), PRef (name2, _) ->
                        let nonterm1 = NonTerminal <| Source.toString name1
                        if not <| nonterminals.Contains nonterm1
                        then
                            nonterminals.Add nonterm1
                        let nonterm2 = NonTerminal <| Source.toString name2
                        if not <| nonterminals.Contains nonterm2
                        then
                            nonterminals.Add nonterm2
                        if not <| crl.ContainsKey (nonterm1, nonterm2)
                        then
                            crl.Add((nonterm1, nonterm2), new ResizeArray<NonTerminal*Probability.T>())
                        if not <| crl.[(nonterm1, nonterm2)].Contains (nonterm, probOne)
                        then
                            crl.[(nonterm1, nonterm2)].Add (nonterm, probOne)                     
                    | _ -> assert false                
                | _ -> assert false

        for key in crl.Keys do
            let list = Seq.toList crl.[key]
            crl_result.Add(key, list)
        for key in srl.Keys do
            let list = Seq.toList srl.[key]
            srl_result.Add(key, list)
        
        let rulesHolder = new RulesHolder(crl_result, srl_result, erl_result)

        recognizeGraph graph squareMatrix rulesHolder nonterminals !S
