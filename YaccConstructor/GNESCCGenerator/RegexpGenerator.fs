﻿//  RegexpGenerator.fs contains functions for list to treee by regexp loading generation
//
//  Copyright 2011 Semen Grigorev <rsdpisuy@gmail.com>
//
//  This file is part of YaccConctructor.
//
//  YaccConstructor is free software:you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

namespace  Yard.Generators.GNESCCGenerator

open Yard.Core.IL
open Yard.Core.IL.Definition
open Yard.Core.IL.Production
open Yard.Core.IL.Rule
open Yard.Generators.GNESCCGenerator.LALR
open Yard.Generators.GNESCCGenerator.FAL
open Yard.Generators.GNESCCGenerator.CommonTypes
open Yard.Generators.GNESCCGenerator.CommonTableFunctions

open Microsoft.FSharp.Text.StructuredFormat
open Microsoft.FSharp.Text.StructuredFormat.LayoutOps

type RegexpGenerator(outPath: string) = 
    class       
        
        let enumerator = new Enumerator()
        
        let ruleToAction = ref []

        let textWriter = TextWriter outPath                                
        let write str = textWriter.Write(str)
         
        let generatePreheader grammarName =
            write "//this file was generated by GNESCC"
            write ("//source grammar:" + grammarName )
            write ("//date:" + System.DateTime.Now.ToString())
            write ""
            write "module GNESCC.Regexp"
            write ""
            write "open Yard.Generators.GNESCCGenerator" 
            write "open System.Text.RegularExpressions"
            write ""
            write "let buildIndexMap kvLst ="    
            write "    let ks = List.map (fun (x:string,y) -> x.Length + 2,y) kvLst"
            write "    List.fold (fun (bl,blst) (l,v) -> bl+l,((bl,v)::blst)) (0,[]) ks"
            write "    |> snd"
            write "    |> dict"
            write ""
            write "let buildStr kvLst ="
            write "    let sep = \";;\""
            write "    List.map fst kvLst "
            write "    |> String.concat sep"
            write "    |> fun s -> \";\" + s + \";\""
            write ""            

        let notMatched expectedType = 
            "| x -> getUnmatched x \"" + expectedType + "\""        

        let rec generateBody (typeToTagMap: System.Collections.Generic.IDictionary<string,int>) groupNum indentSize body =
            let lAltFName = "yardLAltAction"
            let rAltFName = "yardRAltAction"
            let elemFName = "yardElemAction"
            let clsFName = "yardClsAction"
            let optFName = "yardOptAction"
            let indentString l = String.replicate l "    "
                
            match body with 
            | PSeq(elems,expr) ->   

                let genElem i gNum elem =                    
                    generateBody typeToTagMap gNum (indentSize + 1) elem
                    |> fun (body, re, gNum) -> 
                        indentString indentSize + "let e" + i.ToString() + " i =\n"
                        + body
                        ,re
                        ,gNum

                let body,regexp,gns,gn = 
                    List.fold 
                        (fun (bs,rs,gns,gn,i) elem
                            ->
                                let b,re,gn = genElem i gn elem.rule
                                b::bs,re::rs,gn::gns,gn,i+1)
                            ([],[],[],groupNum + 1,0)
                        elems
                    |> fun (b,r,gns,gn,_) ->
                        b |> String.concat "\n"
                        ,
                        r
                        |> List.rev 
                        //|> List.map (fun r -> "(" + r + ")") 
                        |> String.concat ""
                        |> fun r -> "(" + r + ")"
                        ,List.rev gns
                        ,gn
                        
                body + "\n"
                + indentString indentSize + "RESeq [" + (List.mapi (fun i gn -> "e" + string i + " " + gn.ToString()) gns |> String.concat "; ") + "]"
                ,regexp
                ,gn
                 
            | PAlt(alt1,alt2)  ->
                 let genAltItem i gNum elem =
                    generateBody typeToTagMap gNum (indentSize + 1) elem
                    |> fun (body, re, gNum) -> 
                        indentString indentSize + "let e" + i.ToString() + " i =\n"
                        + body
                        ,re
                        ,gNum
                          
                 let lAltB,lAltR,lAltGn = genAltItem 1 (groupNum) alt1
                 let rAltB,rAltR,rAltGn = genAltItem 2 (lAltGn) alt2

                 let regexp = lAltR + "|" + rAltR

                 lAltB + "\n"
                 + rAltB + "\n"
                 + indentString indentSize 
                 + "if elts.[" + (lAltGn |> string) + "].Value = \"\"\n" 
                 + indentString indentSize + "then None, Some (e2"  + rAltGn.ToString() + ")\n"
                 + indentString indentSize + "else Some (e2"  + lAltGn.ToString() + "),None\n"
                 + indentString indentSize + "|> REAlt\n"
                 , regexp
                 , rAltGn
            
            | PRef(x,_) ->
                indentString (indentSize + 1) + "idxValMap.[elts.[" + "i"(*(groupNum + 1).ToString()*) + "].Captures.[0].Index] |> RELeaf"
                , "(;" + typeToTagMap.["NT_" + Source.toString x].ToString() + ";)"
                , groupNum + 1

            | PToken (x) ->
                indentString (indentSize + 1) + "idxValMap.[elts.[" + "i"(*(groupNum + 1).ToString()*) + "].Captures.[0].Index] |> RELeaf"
                , "(;" + typeToTagMap.[ "T_" + Source.toString x].ToString() + ";)"
                , groupNum + 1
               
            | PSome(expr)
            | PMany(expr) as x ->
               let body,re,gn = 
                   generateBody typeToTagMap groupNum (indentSize + 1) expr
               indentString indentSize + "let e i =\n" + body + "\n" 
               + indentString indentSize + "REClosure ["
               + "for c in [0..elts.[" + gn.ToString() + "].Captures.Count-1] -> e c]\n"
               , re + match x with PSome _ -> "+"| PMany _ -> "*" | _ -> failwith "Expected PSome or PMany but get " + x.ToString()
               , gn

//            | POpt(expr) ->
//                 indentString indentSize + "match expr with\n"
//               + indentString indentSize + "| REOpt(opt) -> \n" 
//               + indentString (indentSize + 1) + "let " + optFName + " expr = \n" + (generateBody (indentSize + 2) expr) + "\n"
//               + indentString (indentSize + 1) + "if opt.IsSome then Some (" + optFName + " opt.Value) else None \n"
//               + indentString indentSize + notMatched "REOpt" + "\n"
//               
            | _ -> "NotSupported","NotSupported",groupNum

        let generateRules typeToTagMap rules =            
            let genRule rule =
                let actName =  rule.name
                let arg = "childsLst"
                ruleToAction := (rule.name, actName) :: !ruleToAction
                let body,regexp,gn = generateBody typeToTagMap 0 1 rule.body
                "let " + actName + " " + arg 
                + " = \n    let str = buildStr " + arg 
                + "\n    let idxValMap = buildIndexMap " + arg                
                + "\n    let re = new Regex(\"" + regexp + "\")"
                + "\n    let elts = re.Match(str).Groups" 
                + "\n" + body

            List.map genRule rules

        let generate grammar (typeToTagMap: System.Collections.Generic.IDictionary<string,int>)= 
            generatePreheader grammar.info.fileName            
            generateRules typeToTagMap grammar.grammar |> String.concat "\n" |> write
            List.map 
                (fun r2a -> "(" + (typeToTagMap.["NT_" + fst r2a] |> string) + "," + snd r2a + ")")
                !ruleToAction
            |> String.concat "; "
            |> fun x -> "\nlet ruleToAction = dict [|" + x + "|]\n"
            |> write
            textWriter.CloseOutStream()

        member self.Generate grammar typeToTagMap  = generate grammar typeToTagMap
                
    end