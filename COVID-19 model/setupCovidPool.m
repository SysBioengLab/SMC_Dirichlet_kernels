function pool = setupCovidPool(nWorkers)
% Open a parallel pool for the SMC proposal loop.
%
% Input:   nWorkers  number of workers to request
% Output:  pool      the pool object ([] if none could be opened)
%
% A THREAD pool is tried first (shared memory, so the alive population is not
% copied to each worker every iteration), falling back to a PROCESS pool. With no
% pool, parfor runs serially and only the wall time changes. Reproducibility is
% independent of pool type and worker count: each proposal chunk uses its own
% RandStream substream (see proposeChunk).
% --------------------- Pedro Saa UC 2026 ---------------------------------
if (nargin<1); nWorkers = 4; end

% Clamp to what the local profile allows. Asking for more makes parpool throw
% LocalProfileNumWorkersExceeded, after which parfor silently auto-creates a
% default pool anyway, so the request would be ignored invisibly.
try
    c = parcluster('local');
    if (nWorkers > c.NumWorkers)
        fprintf('Requested %d workers; the local profile allows %d. Using %d.\n',...
                nWorkers,c.NumWorkers,c.NumWorkers);
        nWorkers = c.NumWorkers;
    end
catch
    % no Parallel Computing Toolbox: fall through, parfor will run serially
end

pool = gcp('nocreate');
if ~isempty(pool)
    if (pool.NumWorkers==nWorkers)
        fprintf('Reusing existing %s pool with %d workers.\n',poolKind(pool),pool.NumWorkers);
        return
    end
    delete(pool);                                                           % wrong size, start again
end

try
    pool = parpool('threads',nWorkers);
    fprintf('Opened a THREAD pool with %d workers.\n',pool.NumWorkers);
catch errThread
    if strcmp(errThread.identifier,'parallel:pool:ThreadsNumWorkers')
        % On R2021a a thread pool runs this workload correctly but its size
        % cannot be set; a process pool is used because it honours nWorkers.
        fprintf('Thread pools cannot be limited to %d workers in this release; using a process pool.\n',nWorkers);
    else
        fprintf('Thread pool unavailable (%s). Falling back to a process pool.\n',errThread.identifier);
    end
    try
        pool = parpool('local',nWorkers);
        fprintf('Opened a PROCESS pool with %d workers.\n',pool.NumWorkers);
    catch errProc
        fprintf('No parallel pool could be opened (%s); parfor will run serially.\n',errProc.identifier);
        pool = [];
    end
end
end

%% ------------------------------------------------------------------------
function k = poolKind(pool)
if isa(pool,'parallel.ThreadPool'); k = 'THREAD'; else; k = 'PROCESS'; end
end
