function writeXMLParameterList(filename,mu,beta,lc,delta,nmc)

fp = fopen(filename,'w');

fprintf(fp, '<ParameterList>\n\n');

    fprintf(fp, '\t<ParameterList name="Mesh">\n');
        fprintf(fp, '\t\t<Parameter name="mesh_file" type="string" value="/home/s/staber/Trilinos/nrl/mesh/composite_hexa_32.msh"/>\n');
    fprintf(fp, '\t</ParameterList>\n');

    fprintf(fp, '\t<ParameterList name="nrldata">\n');
        fprintf(fp, '\t\t<Parameter name="pathnrl" type="string" value="/home/s/staber/Trilinos_results/nrl/data/"/>\n');
    fprintf(fp, '\t</ParameterList>\n');

    fprintf(fp, '\t<ParameterList name="Data">\n');
        fprintf(fp, '\t\t<Parameter name="path_to_pts" type="string" value="/home/s/staber/Trilinos_results/nrl/data/xyz.txt"/>\n');
        fprintf(fp, '\t\t<Parameter name="path_to_def" type="string" value="/home/s/staber/Trilinos_results/nrl/data/"/>\n');
    fprintf(fp, '\t</ParameterList>\n');

    fprintf(fp, '\t<ParameterList name="Newton">\n');
        fprintf(fp, '\t\t<Parameter name="delta"             type="double" value="1.0"/>\n');
        fprintf(fp, '\t\t<Parameter name="iterMin"           type="int"    value="2"/>\n');
        fprintf(fp, '\t\t<Parameter name="iterMax"           type="int"    value="10"/>\n');
        fprintf(fp, '\t\t<Parameter name="nbBisMax"          type="int"    value="5"/>\n');
        fprintf(fp, '\t\t<Parameter name="NormFTol"          type="double" value="1e-6"/>\n');
        fprintf(fp, '\t\t<Parameter name="NormFMax"          type="double" value="1e7"/>\n');
        fprintf(fp, '\t\t<Parameter name="eps"               type="double" value="1e-8"/>\n');
        fprintf(fp, '\t\t<Parameter name="success_parameter" type="double" value="2.0"/>\n');
        fprintf(fp, '\t\t<Parameter name="failure_parameter" type="double" value="2.0"/>\n');
        fprintf(fp, '\t\t<Parameter name="number_of_loads"   type="int"    value="1"/>\n');
        fprintf(fp, '\t\t<Parameter name="bc_disp"           type="double" value="1.0"/>\n');
        fprintf(fp, '\t\t<Parameter name="pressure_load"     type="double" value="0.0"/>\n');
        fprintf(fp, '\t\t<Parameter name="tol"               type="double" value="1e-8"/>\n');
    fprintf(fp, '\t</ParameterList>\n\n');

    fprintf(fp, '\t<ParameterList name="Krylov">\n');
        fprintf(fp, '\t\t<Parameter name="solver"          type="string" value="cg"/>\n');
        fprintf(fp, '\t\t<Parameter name="precond"         type="string" value="dom_decomp"/>\n');
        fprintf(fp, '\t\t<Parameter name="subdomain_solve" type="string" value="icc"/>\n');
        fprintf(fp, '\t\t<Parameter name="overlap"         type="int"    value="0"/>\n');
        fprintf(fp, '\t\t<Parameter name="graph_fill"      type="int"    value="0"/>\n');
        fprintf(fp, '\t\t<Parameter name="AZ_tol"          type="double" value="1e-6"/>\n');
        fprintf(fp, '\t\t<Parameter name="AZ_output"       type="int"    value="0"/>\n');
        fprintf(fp, '\t\t<Parameter name="AZ_diagnostics"  type="string" value="all"/>\n');
        fprintf(fp, '\t\t<Parameter name="AZ_reorder"      type="int"    value="1"/>\n');
        fprintf(fp, '\t\t<Parameter name="AZ_conv"         type="string" value="noscaled"/>\n');
    fprintf(fp, '\t</ParameterList>\n\n');

    fprintf(fp, '\t<ParameterList name="TIMooney">\n');
        fprintf(fp, '\t\t<Parameter name="mu1"   type="double" value="%f"/>\n',mu(1));
        fprintf(fp, '\t\t<Parameter name="mu2"   type="double" value="%f"/>\n',mu(2));
        fprintf(fp, '\t\t<Parameter name="mu3"   type="double" value="%f"/>\n',mu(3));
        fprintf(fp, '\t\t<Parameter name="mu4"   type="double" value="%f"/>\n',mu(4));
        fprintf(fp, '\t\t<Parameter name="mu5"   type="double" value="%f"/>\n',mu(5));
        fprintf(fp, '\t\t<Parameter name="beta4" type="double" value="%f"/>\n',beta(1));
        fprintf(fp, '\t\t<Parameter name="beta5" type="double" value="%f"/>\n',beta(2));
    fprintf(fp, '\t</ParameterList>\n\n');

    fprintf(fp, '\t<ParameterList name="Shinozuka">\n');
        fprintf(fp, '\t\t<Parameter name="nmc"    type="int"    value="%d"/>\n',nmc);
        fprintf(fp, '\t\t<Parameter name="order"  type="int"    value="10"/>\n');
        fprintf(fp, '\t\t<Parameter name="lx"     type="double" value="%f"/>\n',lc(1));
        fprintf(fp, '\t\t<Parameter name="ly"     type="double" value="%f"/>\n',lc(2));
        fprintf(fp, '\t\t<Parameter name="delta1" type="double" value="%f"/>\n',delta(1));
        fprintf(fp, '\t\t<Parameter name="delta2" type="double" value="%f"/>\n',delta(2));
        fprintf(fp, '\t\t<Parameter name="delta3" type="double" value="%f"/>\n',delta(3));
        fprintf(fp, '\t\t<Parameter name="delta4" type="double" value="%f"/>\n',delta(4));
    fprintf(fp, '\t</ParameterList>\n\n');

fprintf(fp, '\t</ParameterList>\n');
fclose(fp);

end