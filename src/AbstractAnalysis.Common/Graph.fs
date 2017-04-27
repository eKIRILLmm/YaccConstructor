﻿namespace AbstractAnalysis.Common

open QuickGraph
open System.Runtime.CompilerServices

[<Measure>] type token

[<Measure>] type gssVertex
[<Measure>] type nodeMeasure
[<Measure>] type positionInInput
[<Measure>] type positionInGrammar
[<Measure>] type length
[<Measure>] type leftPosition
[<Measure>] type extension

type LexerEdge<'l ,'br  when 'l: equality> (s, e, t) =
    inherit TaggedEdge<int,Option<'l * 'br>>(s, e, t)
    let l, br =
        match t with
        | Some (l, br) -> Some l, Some br
        | None -> None, None

    member this.BackRef = br
    member this.Label = l

type IParserInput =
    abstract member InitialPositions: array<int<positionInInput>>
        
    abstract member ForAllOutgoingEdges: int<positionInInput> -> (int<token> -> int<positionInInput> -> unit) -> unit
    
    abstract member PositionToString : int -> string

type ParserEdge<'tag>(s, e, t)=
    inherit TaggedEdge<int, 'tag>(s, e, t)
    
type SimpleInputGraph<'tag>(initialVertices : int[], finalVertices : int[], tagToToken : 'tag -> int) = 
    inherit AdjacencyGraph<int, ParserEdge<'tag>>()

    member val InitStates = initialVertices 
    member val FinalStates = finalVertices with get, set
    member val TagToToken = tagToToken with get

    member this.PrintToDot name (tagToString: 'tag -> string) (*(tokenToString : 'token -> string) (numToToken : int -> 'token)*) = 
        use out = new System.IO.StreamWriter (name : string)
        out.WriteLine("digraph AST {")
        out.WriteLine "rankdir=LR"
        for i=0 to this.VertexCount-1 do
            out.Write (i.ToString() + "; ")
        out.WriteLine()
        for i in this.Vertices do
            let edges = this.OutEdges i
            for e in edges do
                let tokenName = e.Tag |> tagToString 
                out.WriteLine (e.Source.ToString() + " -> " + e.Target.ToString() + "[label=\"" + tokenName + "\"]")
        out.WriteLine("}")
        out.Close()      

    new (initial : int, final : int, tagToToken : 'tag -> int) = 
        SimpleInputGraph<_>([|initial|], [|final|], tagToToken)

    new (n : int, tagToToken : 'tag -> int) =
        let allVertices = [|for i in 0 .. n - 1 -> i|]
        SimpleInputGraph<_>(allVertices, allVertices, tagToToken)

    new (initial : int<positionInInput>[], tagToToken : 'tag -> int) = 
        let casted = Array.map(fun x -> int x) initial
        SimpleInputGraph<_>(casted, casted, tagToToken)
 
    interface IParserInput with
        member this.InitialPositions = 
            Array.map(fun x -> x * 1<positionInInput>) this.InitStates

        [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
        member this.ForAllOutgoingEdges curPosInInput pFun =
            let outEdges = int curPosInInput |> this.OutEdges
            outEdges |> Seq.iter
                (fun e -> pFun ((this.TagToToken e.Tag) * 1<token>) (e.Target * 1<positionInInput>))

        member this.PositionToString (pos : int) =
            sprintf "%i" pos


type LinearInput (initialPositions, input:array<int<token>>) =
    interface IParserInput with
        member x.PositionToString(pos: int): string = 
            sprintf "%i" pos

        member this.InitialPositions = initialPositions
        
        [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
        member this.ForAllOutgoingEdges curPosInInput pFun =
            if int curPosInInput < input.Length
            then pFun input.[int curPosInInput] (curPosInInput + 1<positionInInput>)

    member this.Input = input

    new (input:array<int<token>>) = LinearInput ([|0<positionInInput>|], input)

type LinearIputWithErrors(input: int<token> array, errorTag) = 
    interface IParserInput with
        member x.PositionToString(pos: int): string = 
            sprintf "%i" pos

        member this.InitialPositions = [|0<positionInInput>|]
        
        [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
        member this.ForAllOutgoingEdges curPosInInput pFun =
            if int curPosInInput < input.Length
            then 
                pFun input.[int curPosInInput] (curPosInInput + 1<positionInInput>)
                pFun errorTag (curPosInInput + 1<positionInInput>)

    member this.Input = input

type LabelledVertex<'tagType> (pos, tag: 'tagType) =
    member this.Pos = pos
    member this.Tag = tag

type TaggedEdgeWithPos<'vertexType, 'tagType> (pos, f, t, tag: 'tagType) =
    inherit TaggedEdge<'vertexType, 'tagType>(f, t,tag)
    member this.Pos = pos

type EdgeAndVertex<'tagType> = 
    | Edge of TaggedEdgeWithPos<LabelledVertex<'tagType>, 'tagType>
    | Vertex of LabelledVertex<'tagType>

type GraphLabelledVertex<'tagType> (initialVertices : int[], finalVertices : int[], tagToToken : 'tagType -> int) = 
    inherit AdjacencyGraph<LabelledVertex<'tagType>, TaggedEdgeWithPos<LabelledVertex<'tagType>, 'tagType>>()
    let dictEdgeVert = new System.Collections.Generic.Dictionary<int, EdgeAndVertex<'tagType>>()
    member this.AddVerticesAndEdgeRangeWithDict (edges: TaggedEdgeWithPos<LabelledVertex<'tagType>, 'tagType>[]) = 
        this.AddVerticesAndEdgeRange edges |> ignore
        for e in edges do
            dictEdgeVert.Add (e.Pos, Edge(e))
            if not (dictEdgeVert.ContainsKey e.Source.Pos)
            then 
                dictEdgeVert.Add (e.Source.Pos, Vertex(e.Source))
            if not (dictEdgeVert.ContainsKey e.Target.Pos)
            then 
                dictEdgeVert.Add (e.Target.Pos, Vertex(e.Target))
        //printfn "elements in dict: %A" dictEdgeVert.Count


    member val InitStates = initialVertices 
    member val FinalStates = finalVertices with get, set
    member val TagToToken = tagToToken with get

    new (initial : int, final : int, tagToToken : 'tagType -> int) = 
        GraphLabelledVertex<_>([|initial|], [|final|], tagToToken)

    new (n : int, tagToToken : 'tagType -> int) =
        let allPositions = [|for i in 0 .. n - 1 -> i|]
        GraphLabelledVertex<_>(allPositions, allPositions, tagToToken)

    new (initial : int<positionInInput>[], tagToToken : 'tagType -> int) = 
        let casted = Array.map(fun x -> int x) initial
        GraphLabelledVertex<_>(casted, casted, tagToToken)

    interface IParserInput with
        member x.PositionToString(pos: int): string = 
            sprintf "%i" pos
        member this.InitialPositions =
            Array.map(fun x -> x * 1<positionInInput>) this.InitStates
        [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
        member this.ForAllOutgoingEdges curPosInInput pFun =
            match dictEdgeVert.TryGetValue (int curPosInInput) with
            | (true, o) -> 
                match o with
                | Edge e -> pFun ((this.TagToToken e.Target.Tag) * 1<token>) (e.Target.Pos * 1<positionInInput>)
                            //printfn "on edge"
                | Vertex v -> 
                    v |> this.OutEdges
                    |> Seq.iter 
                        (fun edg ->
                            pFun ((this.TagToToken edg.Tag) * 1<token>) (edg.Pos * 1<positionInInput>)
                        )
                    //printfn "on vertex"
            | (_,_) -> printfn "Error"