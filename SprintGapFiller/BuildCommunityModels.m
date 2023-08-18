function [ComModel,removedRxns] = BuildCommunityModels(Folder_path, abbr, Umodel, media, tol)
%%INPUT
%       Folder_path: Matlab cell listing the paths to all the microbial
%                    model structures in .mat format. These models should 
%                    not have exchange reactions in them (EX_) and will be
%                    ignored if it has. Lower (*lb) and upper (*ub) bounds 
%                    can be provided for any reaction in the microbial 
%                    models. Reaction ID (*rxns) and metabolite ID (*mets)
%                    should be in same format as in Umodel. Metabolite IDs 
%                    (*mets) should include the compartment info
%                    (Eg: glc_D(e), pyr(c))
%
%       abbr: Matlab cell listing model abbrevations. All rxns and mets 
%             will have this prefix. Must be same order as in Folder_path
%
%       Umodel: COBRA model structure of Universal model. This model should
%               be a superset of all the reactions and this model should 
%               not have any blocked reactions. All the exchange reactions
%               ID should begin with 'EX_'
%               The following fields are required:
%                   * S - `m x n` Stoichiometric matrix
%                   * b  - `m x 1` change in concentration with time
%                   * c  - `n x 1` Linear objective coefficients
%                   * lb - `n x 1` Lower bounds on net flux
%                   * ub - `n x 1` Upper bounds on net flux
%                   * mets - metabolite IDs
%                   * rxns - reaction IDs
%
%       media: matlab structure with fields
%              *exc_rxns: list of exchange reactions (must be in same 
%                         format as Umodel.rxns)
%              *lb: Lower bounds of corresponding reactions in exc_rxns
%              *ub: Upper bounds of corresponding reactions in exc_rxns
%
%       tol: minimum absolute flux value that every reaction in the 
%            community model should carry (default: 1e-4)

%%OUTPUT
%       ComModel: The consistent community model
%       
%       removedRxns: Reactions that were removed because of inconsistency 
%                    (or) inability to carry the minimum flux (tol) in the given
%                    community and media conditions

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') 
    tol = 1e-4;
end

n_models = numel(Folder_path); % number of models
Umodel_new=Umodel;
exc_rxns=Umodel_new.rxns(startsWith(Umodel_new.rxns,'EX_'));
exc_rxnFormulas = printRxnFormula(Umodel_new,'rxnAbbrList',exc_rxns,'printFlag',false);
Umodel_new=removeRxns(Umodel_new,exc_rxns); % remove all exchange rxns
S=[];lb=[];ub=[];c=[];b=[];rxns=[];mets=[];core=[];
for i=1:n_models
    load(Folder_path{i})
    core = [core;ismember(Umodel_new.rxns,model.rxns)];
    new_rxns = cellfun(@(x)rename_rxns(x,abbr{i}),Umodel_new.rxns,'uni',false);
    rxns = [rxns;new_rxns];
    new_mets=cellfun(@(x)rename_mets(x,abbr{i}),Umodel_new.mets,'uni',false);
    mets=[mets;new_mets];
    S = blkdiag(S,Umodel_new.S);c=[c;Umodel_new.c];b=[b;Umodel_new.b];
    new_lb = Umodel_new.lb;new_ub = Umodel_new.ub;
    [loca,locb] = ismember(model.rxns,Umodel_new.rxns);
    locb = locb(locb~=0);
    new_lb(locb)=model.lb(loca);new_ub(locb)=model.ub(loca);
    lb=[lb;new_lb];ub=[ub;new_ub];
end

% merging the extracellular metabolite rows and removing the extra
% metabolites
[Umets,~,ix] = unique(mets);
counts = accumarray(ix,1).';
counts = counts';
for j=1:numel(counts)
    if counts(j)>1
        ids = find(ismember(mets,Umets{j}));
        S(ids(1),:)=sum(S(ids,:),1);
        S(ids(2:end),:)=[];
        b(ids(2:end))=[];
        mets(ids(2:end))=[];
    end
end
    
ComModel=struct();
ComModel.S=S;ComModel.lb=lb;
ComModel.ub=ub;ComModel.c=c;
ComModel.b=b;ComModel.mets=mets;
ComModel.rxns=rxns;

for i=1:numel(media.exc_rxns)
    id  = find(ismember(exc_rxns,media.exc_rxns{i}));
    ComModel=addReaction(ComModel,exc_rxns{id},'reactionFormula',exc_rxnFormulas{id},...
        'lowerBound',media.lb(i),'upperBound',media.ub(i));
end
core = [core;ones(numel(media.exc_rxns),1)];

a=speedcc(ComModel,tol);

removedRxns = ComModel.rxns(setdiff([1:numel(ComModel.rxns)],a));
community_Umodel = removeRxns(ComModel,removedRxns);
removedRxns = intersect(removedRxns,ComModel.rxns(find(core)));
core =core(a);
ComModel = speedcore(community_Umodel,find(core),tol);

end

function rxns = rename_rxns(a,ABR)
rxns = [ABR,'_',a];
end

function mets = rename_mets(a,ABR)
if ~strcmp(a(end-1),'e')&~contains(a,'biomass')
    mets = [ABR,'_',a];
else
    mets = a;
end
end