function setFigureProperties(varargin)
    % setFigureProperties(NAME,VALUE,...)
    % sets default figure settings I like. Any figure properties can be overrided.

    if isempty(varargin) || ~ishandle(varargin{1})
        fh = gcf;
    else
        fh = varargin{1};
        varargin = varargin(2:end);
    end

    set(fh,'defaultaxesfontsize',14);
    set(fh,'defaultaxesfontname','Arial');
    figSize =[680 228 810 615]; 
    set(fh,'Position',figSize);
    if ~isempty(varargin)
        set(fh,varargin{:});
    end

