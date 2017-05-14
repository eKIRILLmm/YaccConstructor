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

type Pos =
    | Edge = 1
    | Vertex = 0

type GraphLabelledVertex<'tagType when 'tagType : equality> (initialVertices : 'tagType[], finalVertices : 'tagType[], tagToToken : 'tagType -> int) = 
    inherit AdjacencyGraph<'tagType, TaggedEdge<'tagType, 'tagType>>()

    let eMap = new System.Collections.Generic.Dictionary<_,_>()
    let vMap = new System.Collections.Generic.Dictionary<_,_>()
    let vBackMap = new ResizeArray<_>()
    let eBackMap = new ResizeArray<_>()

    let packPosition edge (position: Pos) = 
        if position = Pos.Vertex then ((1 <<< 31) ||| edge) * 1<positionInInput>
        else edge * 1<positionInInput>
    let isVertexOrEdge (position : int<positionInInput>) =
        if int position < 0 then Pos.Vertex
        else Pos.Edge
    let getId (packedValue : int<positionInInput>) = int packedValue &&& 0x7FFFFFFF

    member val InitStates = initialVertices 
    member val FinalStates = finalVertices with get, set
    member val TagToToken = tagToToken with get

    member this.AddEdges (edges: TaggedEdge<'tagType, 'tagType>[]) = 
        this.AddVerticesAndEdgeRange edges |> ignore
        this.Vertices
        |> Seq.iteri (fun i v ->
            vMap.Add(v, i)
            vBackMap.Add v
        )
        let mutable eId = 0
        for e in edges do
            eMap.Add(e, eId) 
            eId <- eId + 1
            eBackMap.Add e


    interface IParserInput with
        member this.InitialPositions = 
            Array.map(fun x -> 
                match (vMap.TryGetValue x) with 
                | (true, v) -> v * 1<positionInInput>
                | (false, v) -> failwithf "There is no vertex %A" v
            ) this.InitStates

        [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
        member this.ForAllOutgoingEdges curPosInInput pFun =
            if (isVertexOrEdge curPosInInput) = Pos.Edge then 
                let e = eBackMap.[getId curPosInInput]
                pFun (this.TagToToken e.Target * 1<token>) (packPosition (vMap.[e.Target]) Pos.Vertex)
            else 
                let outEdges =  vBackMap.[getId curPosInInput] |> this.OutEdges
                outEdges |> Seq.iter
                    (fun e -> pFun ((this.TagToToken e.Tag) * 1<token>) (eMap.[e] * 1<positionInInput>))
//
//        member this.PositionToString (pos : int) =
//            sprintf "%i" pos

         