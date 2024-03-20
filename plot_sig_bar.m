function pl = plot_sig_bar(hdat, mask, vpos, height, Color)
%% plot significance bars for 1D data such as power spectrum in currently open fig
% hdat: hdat to plot on
% mask: row vector of logicals, 1 is significant
% vpos: yaxis location

if nargin == 2 || isempty(vpos)
  vpos = 0.1;
end
ax = gca;
% vpos = ax.YLim(1) + vpos * abs(diff(ax.YLim)); % convert to relative

if nargin == 2
  height = 3;
  Color = [0 0 0];
end

begsample = find(diff( [0 mask 0] ) == 1);
endsample   = find(diff( [0 mask 0] ) == -1) - 1;
if isempty(begsample)
  return % nothing to plot
end
% if begsample == endsample
%     endsample = endsample + 1; % not sure if this is ideal
% end
% pl = plot( hdat([begsample; endsample]), [vpos vpos], 'Linewidth', height, 'Color', Color, 'Marker', 'none', 'Linestyle', '-');

% from ft_plot_vector
%   case 'box'
%     % find the sample number where the highlight begins and endsample
%     begsample = find(diff([0 highlight 0])== 1);
%     endsampleample = find(diff([0 highlight 0])==-1)-1;
%     for i=1:length(begsample)
%       begx = hdat(begsample(i));
%       endx = hdat(endsampleample(i));
%       ft_plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', [.6 .6 .6], 'edgecolor', 'none', 'parent', parent);
%     end

for i=1:length(begsample)
  begx = hdat(begsample(i));
  endx = hdat(endsample(i));
  if begx==endx
    begx = begx - 0.5*mean(diff(hdat));
    endx = endx + 0.5*mean(diff(hdat));    
%     plot(begx, vpos, 'square', 'MarkerFaceColor', Color, 'MarkerEdgeColor', Color)
  end
  pl = plot( [begx; endx], [vpos vpos], 'Linewidth', height, 'Color', Color, 'Marker', 'none', 'Linestyle', '-');
    %         ft_plot_box([begx endx vpos-height/2 vpos+height/2], 'facecolor', Color *0.6, 'edgecolor', 'none');
end



% if isempty(vpos)
%     ax = gca;
%     vpos = ax.YLim(1) + 0.1 * abs(diff(ax.YLim));
% end
%
% begsample = find(diff( [0 mask] ) == 1);
% endsample   = find(diff( [mask 0] ) == -1);
% if isempty(begsample)
%     return % nothing to plot
% end
% if begsample == endsample
%     endsample = endsample + 1; % not sure if this is ideal
% end
% pl = plot( hdat([begsample; endsample]), [vpos vpos], 'height', height, 'Color', Color, 'Marker', 'none', 'Linestyle', '-');
