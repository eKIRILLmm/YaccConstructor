﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using JetBrains.ActionManagement;
using JetBrains.Application.DataContext;
using JetBrains.DataFlow;
using JetBrains.UI.Application;
using JetBrains.UI.ToolWindowManagement;

namespace YC.ReSharper.AbstractAnalysis.Plugin.GraphCodeWindow
{
    [ActionHandler("YC.ReSharper.AbstractAnalysis.Plugin.GraphCodeToolWindow")]
    public class GraphCodeToolWindowAction : IActionHandler
    {
        public bool Update(IDataContext context, ActionPresentation presentation, DelegateUpdate nextUpdate)
        {
            return true;
        }

        public void Execute(IDataContext context, DelegateExecute nextExecute)
        {
            var descriptor = DataConstantsExtensions.GetComponent<GraphCodeToolWindow>(context);
            var manager = DataConstantsExtensions.GetComponent<ToolWindowManager>(context);
            var lifetime = DataConstantsExtensions.GetComponent<Lifetime>(context);
            var uiApplication = DataConstantsExtensions.GetComponent<UIApplication>(context);
            var registrar = new GraphCodeWindowRegistrar(lifetime, manager, descriptor, uiApplication);
            registrar.Show();
        }
    }
}
