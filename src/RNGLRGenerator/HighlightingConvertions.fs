﻿module HighlightingConvertions

open Yard.Core.IL
open Yard.Core.IL.Production
open Yard.Generators.RNGLR

let toClassName (str : string) = 
//    let onlyFirstLetterToUpper (word : string) = 
//        let symbols = [| 
//                        for i = 0 to word.Length - 1 do
//                            if i = 0 
//                            then yield System.Char.ToUpper word.[0]
//                            else yield System.Char.ToLower word.[i] 
//                      |] 
//        new System.String(symbols)
//    
//    str.Split('_')
//    |> Array.map onlyFirstLetterToUpper
//    |> System.String.Concat 
        let symbols = [| 
                        for i = 0 to str.Length - 1 do
                            if i = 0 
                            then yield System.Char.ToUpper str.[0]
                            else yield str.[i] 
                      |] 
        new System.String(symbols)

let litToClassName (lit : string) = 
    toClassName <| lit.ToLower()

let getLeafSemantic leaf isTok = 
    let res = new System.Text.StringBuilder()
    let inline print (x : 'a) =
        Printf.kprintf (fun s -> res.Append s |> ignore) x

    let inline printBr (x : 'a) =
        Printf.kprintf (fun s -> res.Append(s).Append('\n') |> ignore) x

    printBr "let pos = snd <| _rnglr_var_0"
    printBr "let ranges = calculatePos pos"
    printBr "let value = fst <| _rnglr_var_0"

    if isTok
    then printBr "new %sTermNode(\"%s\", value.ToString(), ranges) :> ITreeNode" <| toClassName leaf <| leaf
    else printBr "new %sLitNode(\"%s\", value.ToString(), ranges) :> ITreeNode"  <| litToClassName leaf <| leaf

    res.ToString()
                               
let getNodeSemantic parent children = 
    let res = new System.Text.StringBuilder()
    let inline print (x : 'a) =
        Printf.kprintf (fun s -> res.Append s |> ignore) x

    let inline printBr (x : 'a) =
        Printf.kprintf (fun s -> res.Append(s).Append('\n') |> ignore) x

    let inline printBrInd num (x : 'a) =
        print "%s" (String.replicate (num <<< 2) " ")
        printBr x

    printBrInd 0 "let parent = new %sNonTermNode(\"%s\")" <| toClassName parent <| parent
    printBrInd 0 "let children (*: ITreeNode list*) = %A" children
    printBrInd 0 "addSemantic parent children"
    res.ToString()

let highlightingConvertions (def : Definition.t<Source.t, Source.t>) = 
    let rules = def.grammar.Head.rules    
    let literalToName lit = 
        let indexator = new Indexator(rules, true)
        lit
        |> indexator.literalToIndex
        |> indexator.getLiteralName
        

    let termsAndLitsList = ref []

    let createNewBinding (numOpt : int ref option) = 
        let newBinding = 
            if numOpt.IsSome 
            then 
                incr numOpt.Value
                Some <| new Source.t ("h" + numOpt.Value.Value.ToString())
            else None
        newBinding

    let getNewElem newBinding refRule : elem<Source.t, Source.t> = 
        
        let newElem : elem<Source.t, Source.t> = 
            {
                binding = newBinding
                checker = None
                omit = false
                rule = refRule
            }
        
        newElem

    let changeRule (oldRule : Rule.t<_,_>) (elemList : elem<Source.t, Source.t> list) (bindings : Source.t list) = 
        let actionCode = getNodeSemantic oldRule.name.text bindings
        let newRule : Rule.t<Source.t, Source.t> = 
            {
                name = oldRule.name
                args = []
                body = PSeq(elemList, Some <| new Source.t(actionCode), None)
                isStart = oldRule.isStart
                isPublic = oldRule.isPublic
                metaArgs = []
            }
        newRule

    let rec processElem (oldElem : elem<Source.t, Source.t>) count = 
        let result = ref []
        let newElem = ref oldElem

        let inline createHighlightRefRule text = 
            t.PRef(new Source.t("highlight_" + text), None)
        
        match oldElem.rule with
        | t.PSeq (metaList,_,_) -> 
            for meta in metaList do 
                let someElemList = processElem meta count
                result := !result @ someElemList

        | t.PRef (name, _) as oldElemRule -> 
            let newBinding = createNewBinding <| Some count
            newElem := getNewElem <| newBinding <| t.PRef (name, None)
            result := !result @ [!newElem]

        | t.PToken tok -> 
            if not <| List.exists (fun symbol -> fst symbol = tok.text) !termsAndLitsList
            then termsAndLitsList := (tok.text, true) :: !termsAndLitsList
            
            let newBinding = createNewBinding <| Some count
            newElem := getNewElem newBinding <| createHighlightRefRule tok.text
            result := !result @ [!newElem]
        
        | t.PLiteral lit -> 
            let litName = literalToName lit.text
            if not <| List.exists (fun symbol -> fst symbol = litName) !termsAndLitsList
            then termsAndLitsList := (litName, false) :: !termsAndLitsList
        
            let newBinding = createNewBinding <| Some count
            newElem := getNewElem <| newBinding <| createHighlightRefRule litName
            result := !result @ [!newElem]

        | _ -> failwith "Error in highlighting convertions"
        !result

    let processRule (oldRule : Rule.t<Source.t, Source.t>) = 
        let count = ref 0
        match oldRule.body with 
        | t.PSeq(elemList, _, _) -> 
                            let mutable newElemList = []
                            let mutable bindingsList = []
                            for item in elemList do
                                let newElem = processElem item count
                                newElemList <- newElemList @ newElem

                            for newElem in newElemList do
                                bindingsList <- newElem.binding.Value :: bindingsList

                            changeRule <| oldRule <| newElemList <| List.rev bindingsList
        | t.PRef (_, _) -> oldRule
        | _ -> failwith "Error in highlighting convertions"

    let addHighlightRules() = 
        let mutable res = []
        for tok, isTok in !termsAndLitsList do 
            let actionCode, newElem =
                if isTok 
                then new Source.t(getLeafSemantic tok isTok), getNewElem None <| t.PToken (new Source.t(tok))
                else new Source.t(getLeafSemantic tok isTok), getNewElem None <| t.PLiteral (new Source.t(tok))

            let newRule : Rule.t<Source.t, Source.t> = 
                {
                    name = new Source.t("highlight_" + tok)
                    args = []
                    body = PSeq([newElem], Some <| actionCode, None)
                    isStart = false
                    isPublic = false
                    metaArgs = []
                }
            res <- newRule :: res
        res

    let newRules = List.map processRule rules @ addHighlightRules()

    {def with grammar = [{def.grammar.Head with rules=newRules}]}